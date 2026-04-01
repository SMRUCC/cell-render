Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models
Imports SMRUCC.genomics.Model.OperonMapper

Public Class AnnotationSet

#Region "raw blastp output"
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
#End Region

    Public Property operons As AnnotatedOperon()
    Public Property ec_numbers As Dictionary(Of String, ECNumberAnnotation)

    ''' <summary>
    ''' TFBS data is grouped and index by gene id
    ''' </summary>
    ''' <returns></returns>
    Public Property tfbs_hits As Dictionary(Of String, MotifMatch())
    Public Property transcript_factors As BestHit()
    Public Property membrane_proteins As RankTerm()

End Class
