Imports System.IO
Imports System.Runtime.CompilerServices
Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models
Imports SMRUCC.genomics.Metagenomics
Imports SMRUCC.genomics.SequenceModel

Public Class GenBankProject

    Public Property taxonomy As Taxonomy
    Public Property nt As Dictionary(Of String, Integer)
    Public Property genes As Dictionary(Of String, String)
    Public Property proteins As Dictionary(Of String, String)
    Public Property tss_upstream As Dictionary(Of String, String)
    Public Property gene_table As GeneTable()

    ''' <summary>
    ''' ec number blastp hits
    ''' </summary>
    ''' <returns></returns>
    Public Property enzyme_hits As HitCollection()
    ''' <summary>
    ''' operon gene blastn hits
    ''' </summary>
    ''' <returns></returns>
    Public Property operon_hits As HitCollection()
    Public Property transporter As HitCollection()
    Public Property tf_hits As HitCollection()

    Public Property operons As AnnotatedOperon()
    Public Property ec_numbers As Dictionary(Of String, ECNumberAnnotation)

    ''' <summary>
    ''' TFBS data is grouped and index by gene id
    ''' </summary>
    ''' <returns></returns>
    Public Property tfbs_hits As Dictionary(Of String, MotifMatch())
    Public Property transcript_factors As BestHit()
    Public Property membrane_proteins As RankTerm()

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
