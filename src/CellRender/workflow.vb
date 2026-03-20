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

    <ExportAPI("open_motifdb")>
    Public Function open_motifdb(<RRawVectorArgument> file As Object, Optional env As Environment = Nothing) As Object
        Dim s = SMRUCC.Rsharp.GetFileStream(file, IO.FileAccess.Read, env)

        If s Like GetType(Message) Then
            Return s.TryCast(Of Message)
        End If

        Return MotifDatabase.OpenReadOnly(s.TryCast(Of Stream))
    End Function

    <ExportAPI("motif_search")>
    Public Function motif_search(db As MotifDatabase, <RRawVectorArgument> search_regions As Object,
                                 <RRawVectorArgument(TypeCodes.string)>
                                 Optional family As Object = Nothing,
                                 Optional env As Environment = Nothing) As Object

        Dim seqs As IEnumerable(Of FastaSeq) = pipHelper.GetFastaSeq(search_regions, env)
        Dim familyIds As String() = CLRVector.asCharacter(family)

        If seqs Is Nothing Then
            Return RInternal.debug.stop("invalid fasta sequence source for run TFBS motif site search!", env)
        End If

        Dim motifs As Dictionary(Of String, Probability())

        If familyIds.IsNullOrEmpty Then
            motifs = db.LoadMotifs
        Else
            motifs = familyIds _
                .Distinct _
                .ToDictionary(Function(name) name,
                              Function(name)
                                  Return db _
                                      .LoadFamilyMotifs(name) _
                                      .ToArray
                              End Function)
        End If

        Dim tfbs_hits = db.ScanSites(seqs, n_threads:=1, workflowMode:=True)
    End Function

End Module
