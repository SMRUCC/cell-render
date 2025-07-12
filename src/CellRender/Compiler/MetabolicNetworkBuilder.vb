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
    ReadOnly union_hashcode As UnionHashCode()

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

        ' load data into memory for create cache
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

        union_hashcode = UnionHashCode.LoadUniqueHashCodes(compiler.cad_registry)
        substrate_links = links _
            .GroupBy(Function(a) a.kinetic_id) _
            .ToDictionary(Function(a)
                              Return a.Key.ToString
                          End Function,
                          Function(a)
                              Return a.ToArray
                          End Function)

    End Sub

    ''' <summary>
    ''' get reaction model of a non-enzymatic reaction by its id
    ''' </summary>
    ''' <param name="id"></param>
    ''' <returns></returns>
    Private Function PullReactionNoneEnzymatic(id As String) As Reaction
        Dim r As biocad_registryModel.reaction = cad_registry.reaction.where(field("id") = id).find(Of biocad_registryModel.reaction)

        If r Is Nothing Then
            Return Nothing
        End If

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
            Return Nothing
        End If

        Dim sides = compounds _
            .GroupBy(Function(a) a.side) _
            .ToDictionary(Function(a) a.Key,
                            Function(a)
                                Return a _
                                    .Select(Function(c)
                                                Return New CompoundFactor(c.factor, c.molecule_id.ToString)
                                            End Function) _
                                    .ToArray
                            End Function)

        If Not (sides.ContainsKey("substrate") AndAlso sides.ContainsKey("product")) Then
            Return Nothing
        End If

        Return New Reaction With {
            .ID = r.id,
            .bounds = {1, 1},
            .is_enzymatic = False,
            .name = r.name,
            .substrate = sides!substrate,
            .product = sides!product,
            .note = r.note,
            .ec_number = Nothing
        }
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

    ''' <summary>
    ''' a list of the ec-number that annotated from current genome
    ''' </summary>
    Dim ec_reg As New List(Of String)
    ''' <summary>
    ''' the mapping from gene locus id to the ec_number
    ''' </summary>
    Dim ec_link As New List(Of NamedValue(Of String))

    Private Sub LoadEnzymeLinks()
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
    End Sub

    ''' <summary>
    ''' scan all unique reaction, and build network via rules:
    ''' 
    ''' 1. only add the enzyme reaction which is annotated from the genome
    ''' 2. other non-annotated enzyme reaction is treated as the non-enzymatic reaction
    ''' 3. all non-enzymatic reaction is added into the model
    ''' 4. all reaction is treated as reversiable
    ''' </summary>
    ''' <returns></returns>
    Public Function BuildMetabolicNetwork() As MetabolismStructure
        Dim ec_rxn As Dictionary(Of String, Reaction())
        Dim chemical_rxns As New List(Of Reaction)
        Dim biological_rxns As New List(Of Reaction)

        Call LoadEnzymeLinks()

        Dim ec_numbers As Index(Of String) = ec_reg.Indexing
        Dim ec_generic = ec_numbers.Objects _
            .Select(Function(ec) (ec.Trim("-"c, "."c), ec)) _
            .GroupBy(Function(e) e.Item1) _
            .ToDictionary(Function(e) e.Key,
                          Function(e)
                              Return e.Select(Function(i) i.Item2).Distinct.ToArray
                          End Function)
        Dim enzyme_role = cad_registry.getVocabulary("Enzymatic Catalysis", "Regulation Type")
        Dim scan_id As New Index(Of String)

        For Each hash As UnionHashCode In TqdmWrapper.Wrap(union_hashcode)
            For Each id As String In hash.AsEnumerable
                If id Like scan_id Then
                    Continue For
                Else
                    Call scan_id.Add(id)
                End If

                Dim reaction As Reaction = PullReactionNoneEnzymatic(id)

                If reaction Is Nothing Then
                    Continue For
                End If

                ' check enzyme
                Dim ecs = cad_registry.regulation_graph _
                    .where(field("reaction_id") = id,
                           field("role") = enzyme_role) _
                    .distinct _
                    .project(Of String)("term")

                ecs = ecs.Where(Function(eid)
                                    ' genome contains this ec number exactly
                                    Return eid Like ec_numbers
                                End Function) _
                    .JoinIterates(ec_generic _
                         .AsParallel _
                         .Select(Function(ecg)
                                     ' mapping generic enzyme to annotated enzyme?
                                     If ecs.Any(Function(eid) eid.StartsWith(ecg.Key)) Then
                                         Return ecg.Value
                                     Else
                                         Return {}
                                     End If
                                 End Function) _
                         .IteratesALL) _
                    .Distinct _
                    .ToArray

                If ecs.IsNullOrEmpty Then
                    ' is non-enzymatic reaction
                    Call chemical_rxns.Add(reaction)
                Else
                    ' is enzymatic reaction annotated in current genome
                    reaction.bounds = {10, 10}
                    reaction.is_enzymatic = True
                    reaction.ec_number = ecs
                    reaction.compartment = "Intracellular"

                    Call biological_rxns.Add(reaction)
                End If

                ' hash code is the unique hashcode of the duplicated
                ' reaction models
                ' we has created a reaction object for current hash code
                ' so NO NEEDS for scan other duplicated models
                Exit For
            Next
        Next

        ec_rxn = biological_rxns _
            .Select(Function(a)
                        Return a.ec_number.Select(Function(id) (ec_number:=id, a))
                    End Function) _
            .IteratesALL _
            .GroupBy(Function(r) r.ec_number) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a.Select(Function(g) g.a).ToArray
                          End Function)

        Return New MetabolismStructure With {
            .compounds = PullCompounds(ec_rxn, chemical_rxns.ToArray).ToArray,
            .reactions = New ReactionGroup With {
                .enzymatic = biological_rxns.ToArray,
                .etc = chemical_rxns.ToArray
            },
            .enzymes = queryEnzymes(ec_link, ec_rxn).ToArray
        }
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
                .proteinID = gene.Key,
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
