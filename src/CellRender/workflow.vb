Imports System.IO
Imports CellBuilder
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.SequenceModel.FASTA
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Components
Imports SMRUCC.Rsharp.Runtime.Interop
Imports SMRUCC.Rsharp.Runtime.Vectorization
Imports RInternal = SMRUCC.Rsharp.Runtime.Internal

''' <summary>
''' annotation workflow
''' </summary>
<Package("workflow")>
Module workflow

    <ExportAPI("tss_upstream")>
    <RApiReturn(GetType(FastaSeq))>
    Public Function tss_upstream(proj As GenBankProject) As Object
        Return proj _
            .DumpTSSUpstreamFasta _
            .ToArray
    End Function

End Module
