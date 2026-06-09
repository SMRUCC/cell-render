Imports System.Runtime.CompilerServices
Imports Microsoft.VisualBasic.Linq
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.BLASTOutput.XmlFile
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models

Public Module Annotation

    <Extension>
    Private Iterator Function Parse(enzymes As IEnumerable(Of IGrouping(Of String, Hit)), queryName As String) As IEnumerable(Of ECNumberAnnotation)
        For Each number In enzymes
            Dim sources = number.Select(Function(prot) prot.tag.GetTagValue(" "c, trim:=True)).ToArray
            Dim ec_number As String = number.Key
            Dim total As Double = Aggregate hit As Hit
                                  In number
                                  Into Sum(hit.score * hit.identities * hit.positive)
            Dim source_id As String() = sources.Select(Function(i) i.Name).Distinct.ToArray
            Dim proteinName As IGrouping(Of String, String) = sources _
                .Select(Function(i) i.Value) _
                .GroupBy(Function(i) i) _
                .OrderByDescending(Function(i) i.Count) _
                .First

            Yield New ECNumberAnnotation With {
                .EC = ec_number,
                .Score = total,
                .SourceIDs = source_id,
                .gene_id = queryName,
                .proteinName = proteinName.Key
            }
        Next
    End Function

    <Extension>
    Public Function AssignECNumber(enzymeHits As HitCollection) As ECNumberAnnotation
        Dim enzymes = enzymeHits.AsEnumerable.GroupBy(Function(a) a.hitName.Split("|"c).First).ToArray
        Dim enzyme_scores As ECNumberAnnotation() = enzymes.Parse(enzymeHits.QueryName).ToArray

        If enzyme_scores.Length = 0 Then
            Return Nothing
        ElseIf enzyme_scores.Length = 1 Then
            Return enzyme_scores(0)
        Else
            Return enzyme_scores _
                .OrderByDescending(Function(a) a.Score) _
                .First
        End If
    End Function

    <Extension>
    Private Iterator Function ParseGroups(transcript_factors As IEnumerable(Of IGrouping(Of String, Hit)), queryName As String) As IEnumerable(Of BestHit)
        For Each a In transcript_factors
            Dim ec_number = a.Key
            Dim total As Double = Aggregate hit As Hit In a Into Sum(hit.score * hit.identities * hit.positive)
            Dim source_id As String() = a.Select(Function(prot) prot.hitName.GetTagValue("|"c).Value).Distinct.ToArray

            Yield New BestHit With {
                .HitName = ec_number,
                .hit_length = a.Count,
                .score = total,
                .description = source_id.JoinBy(", "),
                .QueryName = queryName,
                .identities = a.Average(Function(i) i.identities),
                .evalue = a.Average(Function(i) i.evalue),
                .positive = a.Average(Function(i) i.positive)
            }
        Next
    End Function

    <Extension>
    Public Function AssignTFFamilyHit(tfhits As HitCollection) As BestHit
        Dim transcript_factors = tfhits.AsEnumerable.GroupBy(Function(a) a.hitName.Split(" "c).First).ToArray
        Dim tf_scores As BestHit() = transcript_factors.ParseGroups(tfhits.QueryName).ToArray

        If tf_scores.Length = 0 Then
            Return Nothing
        ElseIf tf_scores.Length = 1 Then
            Return tf_scores(0)
        Else
            Return tf_scores _
                .OrderByDescending(Function(a) a.score) _
                .First
        End If
    End Function

End Module

