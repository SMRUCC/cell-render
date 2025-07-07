Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar
Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.ComponentModel.DataSourceModel
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Math.Scripting.MathExpression
Imports Microsoft.VisualBasic.Serialization.JSON
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2

Public Class MetabolicNetworkBuilder

    ReadOnly compiler As Compiler
    ReadOnly chromosome As replicon
    ReadOnly substrate_links As Dictionary(Of String, biocad_registryModel.kinetic_substrate())

    Public ReadOnly Property cad_registry As biocad_registry
        Get
            Return compiler.cad_registry
        End Get
    End Property

    Sub New(compiler As Compiler, chromosome As replicon)
        Me.compiler = compiler
        Me.chromosome = chromosome

        Dim links As New List(Of biocad_registryModel.kinetic_substrate)
        Dim page As biocad_registryModel.kinetic_substrate()
        Dim page_size As Integer = 10000

        For i As Integer = 0 To 100000
            page = compiler.cad_registry.kinetic_substrate _
                .limit(i * page_size, page_size) _
                .select(Of biocad_registryModel.kinetic_substrate)

            If page.IsNullOrEmpty Then
                Exit For
            Else
                Call links.AddRange(page)
            End If
        Next

        substrate_links = links _
            .GroupBy(Function(a) a.kinetic_id) _
            .ToDictionary(Function(a)
                              Return a.Key.ToString
                          End Function,
                          Function(a)
                              Return a.ToArray
                          End Function)
    End Sub

    Private Iterator Function PullReactionNoneEnzymatic() As IEnumerable(Of Reaction)
        Dim page_size = 2000

        For i As Integer = 1 To Integer.MaxValue
            Dim offset As Integer = (i - 1) * page_size
            Dim q = cad_registry.reaction _
                .left_join("regulation_graph") _
                .on(field("`regulation_graph`.reaction_id") = field("`reaction`.id")) _
                .where(field("term").is_nothing) _
                .limit(offset, page_size) _
                .select(Of biocad_registryModel.reaction)("`reaction`.*")

            If q.IsNullOrEmpty Then
                Exit For
            End If

            For Each r As biocad_registryModel.reaction In TqdmWrapper.Wrap(q)
                Dim compounds = cad_registry.reaction_graph _
                    .left_join("vocabulary") _
                    .on(field("`vocabulary`.id") = field("role")) _
                    .where(field("reaction") = r.id) _
                    .select(Of reaction_view)("reaction AS reaction_id",
                                              "molecule_id",
                                              "db_xref",
                                              "term AS side",
                                              "factor")

                If compounds.IsNullOrEmpty OrElse compounds.Any(Function(c) c.molecule_id = 0) Then
                    Continue For
                End If

                Dim sides = compounds _
                    .GroupBy(Function(a) a.side) _
                    .ToDictionary(Function(a) a.Key,
                                    Function(a)
                                        Return a _
                                            .Select(Function(c)
                                                        Return New CompoundFactor(c.factor, c.molecule_id)
                                                    End Function) _
                                            .ToArray
                                    End Function)

                If Not (sides.ContainsKey("substrate") AndAlso sides.ContainsKey("product")) Then
                    Continue For
                End If

                Yield New Reaction With {
                    .ID = r.id,
                    .bounds = {1, 1},
                    .is_enzymatic = False,
                    .name = r.name,
                    .substrate = sides!substrate,
                    .product = sides!product
                }
            Next
        Next
    End Function

    Private Iterator Function FillReactions(ec_rxn As Dictionary(Of String, Reaction())) As IEnumerable(Of Reaction)
        For Each rxn As Reaction In TqdmWrapper.Wrap(ec_rxn.Values _
            .IteratesALL _
            .GroupBy(Function(r) r.ID) _
            .Select(Function(r) r.First) _
            .ToArray)

            Dim reaction = cad_registry.reaction _
                .where(field("id") = rxn.ID) _
                .find(Of biocad_registryModel.reaction)

            rxn.name = reaction.name
            rxn.note = reaction.note

            Yield rxn
        Next
    End Function

    Private Iterator Function PullCompounds(pool As Dictionary(Of String, Reaction()), etc As Reaction()) As IEnumerable(Of Compound)
        Dim all = pool.Values _
            .IteratesALL _
            .Select(Function(rxn) rxn.substrate.JoinIterates(rxn.product)) _
            .IteratesALL _
            .JoinIterates(etc.Select(Function(r) r.substrate.JoinIterates(r.product)).IteratesALL) _
            .GroupBy(Function(c) c.compound) _
            .ToArray
        Dim kegg_id As biocad_registryModel.db_xrefs()

        For Each ref As IGrouping(Of String, CompoundFactor) In TqdmWrapper.Wrap(all)
            Dim mol = cad_registry.molecule _
                .where(field("id") = ref.Key) _
                .find(Of biocad_registryModel.molecule)
            Dim compound As New Compound With {
                .ID = mol.id,
                .name = mol.name,
                .mass0 = 10
            }

            kegg_id = cad_registry.db_xrefs _
                .where(field("db_key") = compiler.kegg_term,
                       field("obj_id") = mol.id) _
                .select(Of biocad_registryModel.db_xrefs)

            If Not kegg_id Is Nothing Then
                compound.kegg_id = kegg_id _
                    .Select(Function(d) d.xref) _
                    .Where(Function(id) id.IsPattern("C\d+")) _
                    .Distinct _
                    .ToArray
            End If

            Yield compound
        Next
    End Function

    Public Function BuildMetabolicNetwork() As MetabolismStructure
        Dim ec_reg As New List(Of String)
        Dim ec_link As New List(Of NamedValue(Of String))
        Dim ec_rxn As Dictionary(Of String, Reaction())
        Dim none_enzymatic = PullReactionNoneEnzymatic().ToArray

        For Each t_unit As TranscriptUnit In TqdmWrapper.Wrap(chromosome.operons)
            For Each gene As gene In t_unit.genes
                ' current gene is not a CDS encoder
                ' skip this rna gene
                If gene.amino_acid Is Nothing Then
                    Continue For
                End If

                Dim ec_numbers As String() = cad_registry.molecule _
                    .left_join("db_xrefs") _
                    .on(field("`db_xrefs`.obj_id") = field("`molecule`.id")) _
                    .where(field("`molecule`.id") = gene.protein_id,
                           field("db_key") = compiler.ec_number) _
                    .distinct() _
                    .project(Of String)("xref")
                Dim prot_id = gene.locus_tag

                Call ec_reg.AddRange(ec_numbers)
                Call ec_link.AddRange(From ec As String
                                      In ec_numbers
                                      Select New NamedValue(Of String)(prot_id, ec))
            Next
        Next

        ec_rxn = loadEnzymeReactions(ec_reg) _
            .ToDictionary(Function(a) a.ec_number,
                          Function(a)
                              Return a.rxns
                          End Function)

        Return New MetabolismStructure With {
            .compounds = PullCompounds(ec_rxn, none_enzymatic).ToArray,
            .reactions = New ReactionGroup With {
                .enzymatic = FillReactions(ec_rxn).ToArray,
                .etc = none_enzymatic
            },
            .enzymes = queryEnzymes(ec_link, ec_rxn).ToArray
        }
    End Function

    Private Iterator Function loadEnzymeReactions(ec_reg As IEnumerable(Of String)) As IEnumerable(Of (ec_number$, rxns As Reaction()))
        For Each ec_number As String In TqdmWrapper.Wrap(ec_reg.Distinct.ToArray)
            Dim ec_generic = ec_number.Trim("-"c, "."c)
            Dim view = cad_registry.regulation_graph _
                .left_join("reaction_graph") _
                .on(field("reaction_graph.reaction") = field("reaction_id")) _
                .left_join("vocabulary") _
                .on(field("vocabulary.id") = field("reaction_graph.role")) _
                .where(field("`regulation_graph`.term") = ec_number Or
                    field("`regulation_graph`.term").instr(ec_generic) = 1) _
                .select(Of reaction_view)("reaction_id",
                                          "molecule_id",
                                          "db_xref",
                                          "vocabulary.term AS side",
                                          "factor")
            Dim reactions As Reaction() = view _
                .Where(Function(c) c.molecule_id > 0 AndAlso Not c.side Is Nothing) _
                .GroupBy(Function(a) a.reaction_id) _
                .Select(Function(rxn)
                            Dim sides = rxn _
                                .GroupBy(Function(a) a.side) _
                                .ToDictionary(Function(a) a.Key,
                                              Function(a)
                                                  Return a _
                                                      .Select(Function(c)
                                                                  Return New CompoundFactor(c.factor, c.molecule_id)
                                                              End Function) _
                                                      .ToArray
                                              End Function)

                            If Not (sides.ContainsKey("substrate") AndAlso sides.ContainsKey("product")) Then
                                Return Nothing
                            End If

                            Return New Reaction With {
                                .ID = rxn.Key,
                                .bounds = {5, 5},
                                .is_enzymatic = True,
                                .name = ec_number,
                                .substrate = sides!substrate,
                                .product = sides!product,
                                .ec_number = {ec_number}
                            }
                        End Function) _
                .Where(Function(r)
                           Return Not r Is Nothing
                       End Function) _
                .ToArray

            Yield (ec_number, reactions)
        Next
    End Function

    Private Iterator Function queryEnzymes(ec_link As IEnumerable(Of NamedValue(Of String)), ec_rxn As Dictionary(Of String, Reaction())) As IEnumerable(Of Enzyme)
        Dim bar As Tqdm.ProgressBar = Nothing

        Call VBDebugger.EchoLine("fetch enzyme and catalysis kinetics data...")

        For Each gene In TqdmWrapper.Wrap(ec_link.GroupBy(Function(a) a.Name).ToArray, bar:=bar)
            Dim ec_str As String() = gene _
                .Select(Function(e) e.Value) _
                .Distinct _
                .ToArray
            Dim rxns = ec_str _
                .Where(Function(id) ec_rxn.ContainsKey(id)) _
                .Select(Function(id) ec_rxn(id)) _
                .IteratesALL _
                .GroupBy(Function(r) r.ID) _
                .ToArray

            Call bar.SetLabel(ec_str.JoinBy(" / "))

            Yield New Enzyme With {
                .geneID = gene.Key,
                .ECNumber = ec_str.JoinBy(" / "),
                .catalysis = rxns _
                    .Select(Function(r) As IEnumerable(Of Catalysis)
                                Return GetKineticsParameters(r)
                            End Function) _
                    .IteratesALL _
                    .ToArray
            }
        Next
    End Function

    Private Iterator Function GetKineticsParameters(r As IGrouping(Of String, Reaction)) As IEnumerable(Of Catalysis)
        ' get ec number for query kinetics law
        Dim ec_id As String() = r.Select(Function(a) a.ec_number).IteratesALL.Distinct.ToArray
        Dim laws = cad_registry.kinetic_law _
            .where(field("ec_number").in(ec_id)) _
            .select(Of biocad_registryModel.kinetic_law)
        ' use substrate network for make confirmed
        Dim hits_any As Boolean = False
        Dim compounds As Index(Of String) = r _
            .Select(Function(a) a.AsEnumerable) _
            .IteratesALL _
            .Select(Function(a) a.compound) _
            .Distinct _
            .Indexing

        For Each law As biocad_registryModel.kinetic_law In laws
            Dim links = substrate_links.TryGetValue(law.id.ToString)

            If links Is Nothing Then
                Continue For
            End If

            For Each meta_link In links
                If meta_link.metabolite_id.ToString Like compounds Then
                    Dim args = law.params.LoadJSON(Of Dictionary(Of String, String))

                    Yield New Catalysis(r.Key) With {
                        .PH = law.pH,
                        .temperature = law.temperature,
                        .parameter = args _
                            .Select(Function(a)
                                        If a.Value.IsNumeric Then
                                            Return New KineticsParameter With {
                                                .name = a.Key,
                                                .value = Val(a.Value)
                                            }
                                        ElseIf a.Value.StartsWith("ENZ") Then
                                            Return New KineticsParameter With {
                                                .name = a.Key,
                                                .isModifier = True,
                                                .target = ec_id.JoinBy("/")
                                            }
                                        Else
                                            Return New KineticsParameter With {
                                                .name = a.Key,
                                                .isModifier = False,
                                                .target = meta_link.metabolite_id
                                            }
                                        End If
                                    End Function) _
                            .ToArray,
                        .formula = New FunctionElement With {.lambda = law.lambda, .name = law.id, .parameters = args.Keys.ToArray}
                    }
                    hits_any = True
                    Exit For
                End If
            Next
        Next

        If Not hits_any Then
            Yield New Catalysis(r.Key) With {
                .PH = 7.0,
                .temperature = 30
            }
        End If
    End Function
End Class
