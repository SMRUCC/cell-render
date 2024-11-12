Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
Imports Microsoft.VisualBasic.ComponentModel.DataSourceModel
Imports Microsoft.VisualBasic.Linq
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular.Vector
Imports SMRUCC.genomics.Metagenomics
Imports [property] = SMRUCC.genomics.GCModeller.CompilerServices.Property

Public Class Compiler

    ReadOnly cad_registry As biocad_registry
    ReadOnly template As GeneTable()
    ReadOnly dna_term As UInteger
    ReadOnly ec_number As UInteger

    Sub New(registry As biocad_registry, genes As GeneTable())
        template = genes
        cad_registry = registry
        dna_term = cad_registry.GetVocabulary("Nucleic Acid").id
        ec_number = cad_registry.GetVocabulary("EC").id
    End Sub

    Private Function BuildGenome() As replicon
        Dim genes As New List(Of TranscriptUnit)

        For Each gene_info As GeneTable In TqdmWrapper.Wrap(template)
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

            If find Is Nothing Then
                Continue For
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

            If Not find_prot Is Nothing Then
                gene.protein_id = find_prot.id
                gene.amino_acid = ProteinComposition _
                    .FromRefSeq(find_prot.sequence, find_prot.xref_id) _
                    .CreateVector
            End If

            Call genes.Add(New TranscriptUnit With {.id = find.id, .genes = {gene}})
        Next

        Return New replicon With {
            .genomeName = "",
            .isPlasmid = False,
            .operons = genes.ToArray,
            .RNAs = {}
        }
    End Function

    Private Function BuildMetabolicNetwork(chromosome As replicon) As MetabolismStructure
        Dim ec_reg As New List(Of String)
        Dim ec_link As New List(Of NamedValue(Of String))
        Dim ec_rxn As New Dictionary(Of String, Reaction())

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
                Dim prot_id = gene.amino_acid.name

                Call ec_reg.AddRange(ec_numbers)
                Call ec_link.AddRange(From ec As String
                                      In ec_numbers
                                      Select New NamedValue(Of String)(prot_id, ec))
            Next
        Next

        For Each ec_number As String In TqdmWrapper.Wrap(ec_reg.Distinct.ToArray)
            Dim view = cad_registry.regulation_graph _
                .left_join("reaction_graph") _
                .on(field("reaction_graph.reaction") = field("reaction_id")) _
                .left_join("vocabulary") _
                .on(field("vocabulary.id") = field("reaction_graph.role")) _
                .where(field("regulation_graph.term") = ec_number) _
                .select(Of reaction_view)("reaction_id", "molecule_id", "db_xref", "vocabulary.term AS side", "factor")
            Dim reactions As Reaction() = view _
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

                            Return New Reaction With {
                                .ID = rxn.Key,
                                .bounds = {1, 1},
                                .is_enzymatic = True,
                                .name = ec_number,
                                .substrate = sides!substrate,
                                .product = sides!product
                            }
                        End Function) _
                .ToArray

            Call ec_rxn.Add(ec_number, reactions)
        Next

        Return New MetabolismStructure With {
            .compounds = PullCompounds(ec_rxn).ToArray,
            .reactions = New ReactionGroup With {
                .enzymatic = ec_rxn.Values _
                    .IteratesALL _
                    .GroupBy(Function(r) r.ID) _
                    .Select(Function(r) r.First) _
                    .ToArray
            },
            .enzymes =
        }
    End Function

    Private Iterator Function PullCompounds(pool As Dictionary(Of String, Reaction())) As IEnumerable(Of Compound)
        Dim all = pool.Values _
            .IteratesALL _
            .Select(Function(rxn) rxn.substrate.JoinIterates(rxn.product)) _
            .IteratesALL _
            .GroupBy(Function(c) c.compound) _
            .ToArray

        For Each ref In TqdmWrapper.Wrap(all)
            Dim mol = cad_registry.molecule.where(field("id") = ref.Key).find(Of biocad_registryModel.molecule)
            Dim compound As New Compound With {
                .ID = mol.id,
                .name = mol.name,
                .mass0 = 10
            }

            Yield compound
        Next
    End Function

    Public Function CreateModel() As VirtualCell
        Dim chromosome As replicon = BuildGenome()
        Dim metabolic As MetabolismStructure = BuildMetabolicNetwork(chromosome)

        Return New VirtualCell With {
            .properties = New [property],
            .taxonomy = New Taxonomy,
            .genome = New Genome With {
                .replicons = {chromosome}
            },
            .metabolismStructure = metabolic
        }
    End Function
End Class
