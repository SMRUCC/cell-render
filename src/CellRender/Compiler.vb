﻿Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
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

    Sub New(registry As biocad_registry, genes As GeneTable())
        template = genes
        cad_registry = registry
        dna_term = cad_registry.GetVocabulary("Nucleic Acid").id
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
                .find(Of gene_molecule)("molecule.xref_id", "sequence")
            Dim gene As New gene(gene_info.Location) With {
                .locus_tag = gene_info.locus_id,
                .product = find.note,
                .nucleotide_base = rna
            }

            If Not find_prot Is Nothing Then
                gene.protein_id = find_prot.xref_id
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

    Public Function CreateModel() As VirtualCell
        Dim chromosome As replicon = BuildGenome()

        Return New VirtualCell With {
            .properties = New [property],
            .taxonomy = New Taxonomy,
            .genome = New Genome With {
                .replicons = {chromosome}
            }
        }
    End Function
End Class
