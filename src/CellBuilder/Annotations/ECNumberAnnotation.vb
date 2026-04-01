Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models

Public Class ECNumberAnnotation

    Public Property gene_id As String
    Public Property EC As String
    Public Property Score As Double
    Public Property SourceIDs As String()

    Public Shared Iterator Function MakeEnzymeTerms(blastp As IEnumerable(Of HitCollection)) As IEnumerable(Of ECNumberAnnotation)
        For Each protein As HitCollection In blastp
            Dim ec As ECNumberAnnotation = protein.AssignECNumber()

            If Not ec Is Nothing Then
                Yield ec
            End If
        Next
    End Function

End Class