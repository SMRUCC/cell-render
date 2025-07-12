Imports System.Runtime.InteropServices
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Assembly.NCBI.GenBank
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.Metagenomics
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
                                    Optional taxid As String = Nothing,
                                    Optional env As Environment = Nothing) As Object

        Dim template As GeneTable()
        Dim taxinfo As Taxonomy = Nothing
        Dim cellular_id As String = "intracellular"

        If genes Is Nothing Then
            Return RInternal.debug.stop("the required template source should not be nothing!", env)
        End If

        If TypeOf genes Is GBFF.File Then
            template = DirectCast(genes, GBFF.File) _
                .EnumerateGeneFeatures(ORF:=False) _
                .AsParallel _
                .Select(Function(gene) gene.DumpExportFeature) _
                .ToArray
            taxinfo = DirectCast(genes, GBFF.File).Source.GetTaxonomy
            cellular_id = If(DirectCast(genes, GBFF.File).Taxon, cellular_id)
        Else
            Dim pull As pipeline = pipeline.TryCreatePipeline(Of GeneTable)(genes, env)

            If pull.isError Then
                Return Message.InCompatibleType(GetType(GBFF.File), genes.GetType, env)
            Else
                template = pull.populates(Of GeneTable)(env).ToArray
            End If
        End If

        If taxid.IsPattern("\d+") Then
            template = template _
               .Select(Function(gene)
                           gene.locus_id = taxid & ":" & gene.locus_id
                           Return gene
                       End Function) _
               .ToArray
        End If

        Using compiler As New Compiler(register, template, CLng(Val(taxid.Match("\d+"))).ToString, cellular_id)
            Dim vcell As VirtualCell = compiler.Compile($"/compile --log ""{logfile}""")

            vcell.taxonomy = taxinfo
            vcell.cellular_id = cellular_id

            Return vcell
        End Using
    End Function

End Module
