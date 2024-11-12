Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.Metagenomics

Public Class Compiler

    ReadOnly cad_registry As biocad_registry
    ReadOnly template As GeneTable()

    Sub New(registry As biocad_registry, genes As GeneTable())
        template = genes
        cad_registry = registry
    End Sub

    Private Function BuildGenome() As replicon

    End Function

    Public Function CreateModel() As VirtualCell
        Dim chromosome As replicon = BuildGenome()

        Return New VirtualCell With {
            .properties = New [Property],
            .taxonomy = New Taxonomy,
            .genome = New Genome With {
                .replicons = {chromosome}
            }
        }
    End Function
End Class
