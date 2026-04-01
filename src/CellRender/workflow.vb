Imports System.IO
Imports CellBuilder
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models
Imports SMRUCC.genomics.SequenceModel.FASTA
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Components
Imports SMRUCC.Rsharp.Runtime.Internal.[Object]
Imports SMRUCC.Rsharp.Runtime.Interop
Imports SMRUCC.Rsharp.Runtime.Vectorization
Imports RInternal = SMRUCC.Rsharp.Runtime.Internal

''' <summary>
''' annotation workflow
''' </summary>
<Package("workflow")>
Module workflow

    ''' <summary>
    ''' extract of the tss upstream location site sequence data
    ''' </summary>
    ''' <param name="proj"></param>
    ''' <returns></returns>
    <ExportAPI("tss_upstream")>
    <RApiReturn(GetType(FastaSeq))>
    Public Function tss_upstream(proj As GenBankProject) As Object
        Return proj _
            .DumpTSSUpstreamFasta _
            .ToArray
    End Function

    ''' <summary>
    ''' extract of the protein fasta sequence data to file
    ''' </summary>
    ''' <param name="proj"></param>
    ''' <param name="file"></param>
    ''' <param name="env"></param>
    ''' <returns></returns>
    <ExportAPI("save_proteins")>
    Public Function save_proteins(proj As GenBankProject, file As Object, Optional env As Environment = Nothing) As Object
        Dim s = SMRUCC.Rsharp.GetFileStream(file, FileAccess.Write, env)

        If s Like GetType(Message) Then
            Return s.TryCast(Of Message)
        End If

        Call proj.DumpProteinFasta(s.TryCast(Of Stream))

        Return Nothing
    End Function

    <ExportAPI("set_blastp_result")>
    Public Function set_blastp_result(<RRawVectorArgument> blastp_hits As Object, proj As GenBankProject, group As String, Optional env As Environment = Nothing) As Object
        Dim pull As pipeline = pipeline.TryCreatePipeline(Of HitCollection)(blastp_hits, env)

        If pull.isError Then
            Return pull.getError
        ElseIf proj.annotations Is Nothing Then
            proj.annotations = New AnnotationSet
        End If

        Dim protein_hits As HitCollection() = pull.populates(Of HitCollection)(env).ToArray

        Select Case Strings.LCase(group)
            Case "ec_number" : proj.annotations.enzyme_hits = protein_hits
            Case "subcellular_location" : proj.annotations.transporter = protein_hits
            Case "transcript_factor" : proj.annotations.tf_hits = protein_hits
            Case Else
                Return RInternal.debug.stop($"unknown annotation group of '{group}'", env)
        End Select

        Return Nothing
    End Function

    <ExportAPI("set_tfbs")>
    Public Function set_tfbs(proj As GenBankProject, <RRawVectorArgument> tfbs As Object, Optional env As Environment = Nothing) As Object
        Dim pull As pipeline = pipeline.TryCreatePipeline(Of MotifMatch)(tfbs, env)

        If pull.isError Then
            Return pull.getError
        ElseIf proj.annotations Is Nothing Then
            proj.annotations = New AnnotationSet
        End If

        proj.annotations.tfbs_hits = pull _
            .populates(Of MotifMatch)(env) _
            .GroupBy(Function(m) m.title) _
            .ToDictionary(Function(m) m.Key,
                          Function(m)
                              Return m.ToArray
                          End Function)

        Return proj
    End Function
End Module
