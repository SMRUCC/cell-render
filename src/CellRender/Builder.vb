Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Assembly.NCBI.GenBank
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Components
Imports SMRUCC.Rsharp.Runtime.Internal.[Object]
Imports SMRUCC.Rsharp.Runtime.Interop
Imports RInternal = SMRUCC.Rsharp.Runtime.Internal

''' <summary>
''' Helper functions for build virtualcell model
''' </summary>
<Package("Builder")>
<RTypeExport("cad_registry", GetType(biocad_registry))>
Public Module Builder

    <ExportAPI("create_modelfile")>
    <RApiReturn(GetType(VirtualCell))>
    Public Function CreateModelFile(register As biocad_registry,
                                    <RRawVectorArgument>
                                    genes As Object,
                                    Optional logfile As String = "./model_compile.log",
                                    Optional env As Environment = Nothing) As Object

        Dim template As GeneTable()

        If genes Is Nothing Then
            Return RInternal.debug.stop("the required template source should not be nothing!", env)
        End If

        If TypeOf genes Is GBFF.File Then
            template = DirectCast(genes, GBFF.File) _
                .EnumerateGeneFeatures(ORF:=False) _
                .AsParallel _
                .Select(Function(gene) gene.DumpExportFeature) _
                .ToArray
        Else
            Dim pull As pipeline = pipeline.TryCreatePipeline(Of GeneTable)(genes, env)

            If pull.isError Then
                Return Message.InCompatibleType(GetType(GBFF.File), genes.GetType, env)
            Else
                template = pull.populates(Of GeneTable)(env).ToArray
            End If
        End If

        Using compiler As New Compiler(register, template)
            Return compiler.Compile($"/compile --log ""{logfile}""")
        End Using
    End Function

End Module
