Imports Microsoft.VisualBasic.CommandLine
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.Metagenomics
Imports [property] = SMRUCC.genomics.GCModeller.CompilerServices.Property

''' <summary>
''' Compiler for the gcmodeller virtual cell model based on the ncbi genbank gene models.
''' </summary>
Public Class Compiler : Inherits Compiler(Of VirtualCell)

    Friend ReadOnly cad_registry As biocad_registry
    Friend ReadOnly template As GeneTable()
    Friend ReadOnly dna_term As UInteger
    Friend ReadOnly ec_number As UInteger
    Friend ReadOnly kegg_term As UInteger
    Friend ReadOnly polypeptide_term As UInteger
    Friend ReadOnly tax_id As String

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="registry"></param>
    ''' <param name="genes">
    ''' A gene set that loaded from the NCBI genbank file.
    ''' </param>
    ''' <remarks>
    ''' The gene set should be a <see cref="GeneTable"/> array that contains all the genes
    ''' information of a genome, such as locus_tag, product, location, etc.
    ''' </remarks>
    Sub New(registry As biocad_registry, genes As GeneTable(), taxid As String)
        Dim terms As BioCadVocabulary = registry.vocabulary_terms

        template = genes.GroupBy(Function(g) g.locus_id) _
            .Select(Function(g) g.First) _
            .ToArray
        cad_registry = registry
        dna_term = terms.gene_term
        ec_number = terms.ecnumber_term
        kegg_term = terms.kegg_term
        polypeptide_term = terms.protein_term
        tax_id = taxid
    End Sub

    Protected Overrides Function CompileImpl(args As CommandLine) As Integer
        Dim chromosome As replicon = New ReplicateBuilder(Me).BuildGenome()
        Dim metabolic As MetabolismStructure = New MetabolicNetworkBuilder(Me, chromosome).BuildMetabolicNetwork()

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
