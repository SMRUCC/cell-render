Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.Metagenomics
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Interop
Imports SMRUCC.Rsharp.Runtime.Vectorization

''' <summary>
''' Helper functions for build virtualcell model
''' </summary>
<Package("Builder")>
Public Module Builder

    <ExportAPI("create_modelfile")>
    <RApiReturn(GetType(VirtualCell))>
    Public Function CreateModelFile(register As biocad_registry.biocad_registry,
                                    <RRawVectorArgument>
                                    genes As Object,
                                    Optional env As Environment = Nothing) As Object

        Dim geneIds As String() = CLRVector.asCharacter(genes)
        Dim chromosome As New replicon

        Return New VirtualCell With {
            .properties = New [Property],
            .taxonomy = New Taxonomy,
            .genome = New Genome With {
                .replicons = {chromosome}
            }
        }
    End Function

End Module
