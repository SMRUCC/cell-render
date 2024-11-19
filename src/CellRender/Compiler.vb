Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar
Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
Imports Microsoft.VisualBasic.CommandLine
Imports Microsoft.VisualBasic.ComponentModel.DataSourceModel
Imports Microsoft.VisualBasic.Linq
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular.Vector
Imports SMRUCC.genomics.Metagenomics
Imports [property] = SMRUCC.genomics.GCModeller.CompilerServices.Property

Public Class Compiler : Inherits Compiler(Of VirtualCell)

    ReadOnly cad_registry As biocad_registry
    ReadOnly template As GeneTable()
    ReadOnly dna_term As UInteger
    ReadOnly ec_number As UInteger
    ReadOnly kegg_term As UInteger
    ReadOnly polypeptide_term As UInteger

    Sub New(registry As biocad_registry, genes As GeneTable())
        template = genes
        cad_registry = registry
        dna_term = cad_registry.GetVocabulary("Nucleic Acid").id
        ec_number = cad_registry.GetVocabulary("EC").id
        kegg_term = cad_registry.GetVocabulary("KEGG").id
        polypeptide_term = cad_registry.GetVocabulary("Polypeptide").id
    End Sub

    Private Function BuildGenome() As replicon
        Dim genes As New List(Of TranscriptUnit)
        Dim rnas As New List(Of RNA)
        Dim bar As Tqdm.ProgressBar = Nothing

        Call VBDebugger.EchoLine("compile of the genome model, pull gene and proteins.")

        ' contains CDS/tRNA/rRNA
        For Each gene_info As GeneTable In TqdmWrapper.Wrap(template, bar:=bar, useColor:=True)
            ' fetch gene information from database
            Dim find As gene_molecule = cad_registry.molecule _
                .left_join("db_xrefs") _
                .on(field("`db_xrefs`.obj_id") = field("`molecule`.id")) _
                .left_join("sequence_graph") _
                .on(field("`sequence_graph`.molecule_id") = field("`molecule`.id")) _
                .where(field("`molecule`.type") = dna_term,
                       field("xref_id") = gene_info.locus_id Or
                       field("xref").in({
                           gene_info.locus_id,
                           gene_info.ProteinId,
                           gene_info.UniprotSwissProt,
                           gene_info.UniprotTrEMBL}, nullFilter:=True)) _
                .find(Of gene_molecule)("`molecule`.id", "xref_id", "name", "note", "sequence")

            ' missing current gene item inside database
            If find Is Nothing Then
                Continue For
            Else
                Call bar.SetLabel(gene_info.ToString)
            End If

            Dim rna = RNAComposition _
                .FromNtSequence(find.sequence, gene_info.locus_id) _
                .CreateVector
            Dim find_prot = cad_registry.molecule _
                .left_join("sequence_graph") _
                .on(field("`sequence_graph`.molecule_id") = field("`molecule`.id")) _
                .where(field("`molecule`.parent") = find.id) _
                .find(Of gene_molecule)("`molecule`.id", "molecule.xref_id", "sequence")
            Dim gene As New gene(gene_info.Location) With {
                .locus_tag = gene_info.locus_id,
                .product = find.note,
                .nucleotide_base = rna
            }

            If find_prot Is Nothing AndAlso Not gene_info.ProteinId.StringEmpty(, True) Then
                find_prot = cad_registry.db_xrefs _
                    .left_join("sequence_graph") _
                    .on(field("`sequence_graph`.molecule_id") = field("obj_id")) _
                    .where(field("type") = polypeptide_term, field("xref") = gene_info.ProteinId) _
                    .find(Of gene_molecule)("molecule_id AS id", "xref AS xref_id", "sequence")
            End If

            If Not find_prot Is Nothing Then
                ' find a protein sequnece
                ' is CDS/ORF
                gene.protein_id = find_prot.id
                gene.amino_acid = ProteinComposition _
                    .FromRefSeq(find_prot.sequence, find_prot.xref_id) _
                    .CreateVector
            Else
                ' no protein sequence could be found
                ' is rRNA or tRNA or other kind of RNA
                Call rnas.Add(New RNA With {
                    .gene = gene.locus_tag,
                    .type = RNAType(gene_info.type),
                    .val = gene_info.commonName
                })
            End If

            Call genes.Add(New TranscriptUnit With {.id = find.id, .genes = {gene}})
        Next

        Return New replicon With {
            .genomeName = "",
            .isPlasmid = False,
            .operons = genes.ToArray,
            .RNAs = rnas.ToArray
        }
    End Function

    Private Shared Function RNAType(s As String) As RNATypes
        Select Case Strings.Trim(s).ToLower
            ' broken data!
            Case "cds" : Return RNATypes.micsRNA
            Case "trna"
                Return RNATypes.tRNA
            Case "rrna"
                Return RNATypes.ribosomalRNA
            Case Else
                Return RNATypes.micsRNA
        End Select
    End Function

    Private Function BuildMetabolicNetwork(chromosome As replicon) As MetabolismStructure
        Dim ec_reg As New List(Of String)
        Dim ec_link As New List(Of NamedValue(Of String))
        Dim ec_rxn As New Dictionary(Of String, Reaction())
        Dim none_enzymatic = PullReactionNoneEnzyme().ToArray

        For Each t_unit As TranscriptUnit In TqdmWrapper.Wrap(chromosome.operons)
            For Each gene As gene In t_unit.genes
                If gene.amino_acid Is Nothing Then
                    Continue For
                End If

                Dim ec_numbers As String() = cad_registry.molecule _
                    .left_join("db_xrefs") _
                    .on(field("`db_xrefs`.obj_id") = field("`molecule`.id")) _
                    .where(field("`molecule`.id") = gene.protein_id,
                           field("db_key") = ec_number) _
                    .distinct() _
                    .project(Of String)("xref")
                Dim prot_id = gene.locus_tag

                Call ec_reg.AddRange(ec_numbers)
                Call ec_link.AddRange(From ec As String
                                      In ec_numbers
                                      Select New NamedValue(Of String)(prot_id, ec))
            Next
        Next

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
                                .product = sides!product
                            }
                        End Function) _
                .Where(Function(r) Not r Is Nothing) _
                .ToArray

            Call ec_rxn.Add(ec_number, reactions)
        Next

        Dim enzymes As New List(Of Enzyme)

        For Each gene In ec_link.GroupBy(Function(a) a.Name)
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

            Call enzymes.Add(New Enzyme With {
                .geneID = gene.Key,
                .ECNumber = ec_str.JoinBy(" / "),
                .catalysis = rxns _
                    .Select(Function(r)
                                Return New Catalysis(r.Key) With {
                                    .PH = 7.0,
                                    .temperature = 30
                                }
                            End Function) _
                    .ToArray
            })
        Next

        Return New MetabolismStructure With {
            .compounds = PullCompounds(ec_rxn, none_enzymatic).ToArray,
            .reactions = New ReactionGroup With {
                .enzymatic = FillReactions(ec_rxn).ToArray,
                .etc = none_enzymatic
            },
            .enzymes = enzymes.ToArray
        }
    End Function

    Private Iterator Function PullReactionNoneEnzyme() As IEnumerable(Of Reaction)
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
                .where(field("db_key") = kegg_term,
                       field("obj_id") = mol.id) _
                .select(Of biocad_registryModel.db_xrefs)

            If Not kegg_id Is Nothing Then
                compound.kegg_id = kegg_id _
                    .Select(Function(d) d.xref) _
                    .Distinct _
                    .ToArray
            End If

            Yield compound
        Next
    End Function

    Protected Overrides Function CompileImpl(args As CommandLine) As Integer
        Dim chromosome As replicon = BuildGenome()
        Dim metabolic As MetabolismStructure = BuildMetabolicNetwork(chromosome)

        m_compiledModel = New VirtualCell With {
            .properties = New [property],
            .taxonomy = New Taxonomy,
            .genome = New Genome With {
                .replicons = {chromosome}
            },
            .metabolismStructure = metabolic
        }

        Return 0
    End Function
End Class
