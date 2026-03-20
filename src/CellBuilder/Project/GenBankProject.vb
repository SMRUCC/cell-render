Imports System.IO
Imports System.Runtime.CompilerServices
Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
Imports Microsoft.VisualBasic.Linq
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.Metagenomics
Imports SMRUCC.genomics.SequenceModel
Imports SMRUCC.genomics.SequenceModel.FASTA

Public Class GenBankProject

    Public Property taxonomy As Taxonomy
    Public Property nt As Dictionary(Of String, Integer)
    Public Property genes As Dictionary(Of String, String)
    Public Property proteins As Dictionary(Of String, String)
    Public Property tss_upstream As Dictionary(Of String, String)
    Public Property gene_table As GeneTable()
    Public Property annotations As AnnotationSet

    <MethodImpl(MethodImplOptions.AggressiveInlining)>
    Public Sub DumpProteinFasta(s As Stream)
        Call FASTA.StreamWriter.WriteList(proteins, s)
    End Sub

    <MethodImpl(MethodImplOptions.AggressiveInlining)>
    Public Sub DumpGeneFasta(s As Stream)
        Call FASTA.StreamWriter.WriteList(genes, s)
    End Sub

    <MethodImpl(MethodImplOptions.AggressiveInlining)>
    Public Sub DumpTSSUpstreamFasta(s As Stream)
        Call FASTA.StreamWriter.WriteList(tss_upstream, s)
    End Sub

    Public Iterator Function DumpTSSUpstreamFasta() As IEnumerable(Of FastaSeq)
        For Each seq In tss_upstream.SafeQuery
            Yield New FastaSeq(seq.Value, title:=seq.Key)
        Next
    End Function

    Public Function ComputeHashCode() As String
        Dim hashcode As String = ""

        hashcode = SequenceHashcode(genes) & SequenceHashcode(proteins)
        hashcode = hashcode.MD5

        Return hashcode
    End Function

    Private Shared Function SequenceHashcode(seqSet As Dictionary(Of String, String)) As String
        Dim sort = seqSet.OrderBy(Function(s) s.Key).ToArray
        Dim hashcode As String = ""

        For Each term In TqdmWrapper.Wrap(sort)
            hashcode = (hashcode & term.Key & term.Value).MD5
        Next

        Return hashcode
    End Function

End Class
