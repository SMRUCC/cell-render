Imports System.Runtime.CompilerServices
Imports Microsoft.VisualBasic.Linq
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models

Public Module Annotation

    <Extension>
    Public Function AssignECNumber(enzymeHits As HitCollection) As ECNumberAnnotation
        Dim enzymes = enzymeHits.AsEnumerable.GroupBy(Function(a) a.hitName.Split("|"c).First).ToArray
        Dim enzyme_scores As ECNumberAnnotation() = enzymes _
            .Select(Function(a)
                        Dim ec_number = a.Key
                        Dim total As Double = Aggregate hit As Hit In a Into Sum(hit.score * hit.identities * hit.positive)
                        Dim source_id As String() = a.Select(Function(prot) prot.hitName.GetTagValue("|"c).Value).Distinct.ToArray

                        Return New ECNumberAnnotation With {
                            .EC = ec_number,
                            .Score = total,
                            .SourceIDs = source_id,
                            .gene_id = enzymeHits.QueryName
                        }
                    End Function) _
            .ToArray

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
    Public Function AssignTFFamilyHit(tfhits As HitCollection) As BestHit
        Dim enzymes = tfhits.AsEnumerable.GroupBy(Function(a) a.hitName.Split("|"c).First).ToArray
        Dim enzyme_scores As BestHit() = enzymes _
            .Select(Function(a)
                        Dim ec_number = a.Key
                        Dim total As Double = Aggregate hit As Hit In a Into Sum(hit.score * hit.identities * hit.positive)
                        Dim source_id As String() = a.Select(Function(prot) prot.hitName.GetTagValue("|"c).Value).Distinct.ToArray

                        Return New BestHit With {
                            .HitName = ec_number,
                            .hit_length = a.Count,
                            .score = total,
                            .description = source_id.JoinBy(", "),
                            .QueryName = tfhits.QueryName,
                            .identities = a.Average(Function(i) i.identities),
                            .evalue = a.Average(Function(i) i.evalue),
                            .positive = a.Average(Function(i) i.positive)
                        }
                    End Function) _
            .ToArray

        If enzyme_scores.Length = 0 Then
            Return Nothing
        ElseIf enzyme_scores.Length = 1 Then
            Return enzyme_scores(0)
        Else
            Return enzyme_scores _
                .OrderByDescending(Function(a) a.score) _
                .First
        End If
    End Function

End Module

