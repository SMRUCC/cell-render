Imports System.Runtime.CompilerServices
Imports CellBuilder
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Linq
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.Interops.NCBI.Extensions
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models

Public Module SetBlastResult

    <Extension>
    Public Sub SetECNumber(proj As GenBankProject, protein_hits As HitCollection())
        proj.annotations.enzyme_hits = protein_hits
        proj.annotations.ec_numbers = ECNumberAnnotation _
            .MakeEnzymeTerms(proj.annotations.enzyme_hits) _
            .ToDictionary(Function(a)
                              Return a.gene_id
                          End Function)
        ' 20260506
        ' update ec numbers of the gene annotation table
        For Each gene As GeneTable In proj.gene_table
            If proj.annotations.ec_numbers.ContainsKey(gene.locus_id) Then
                gene.EC_Number = {proj.annotations.ec_numbers(gene.locus_id).EC}
                gene.COG = proj.annotations.ec_numbers(gene.locus_id).proteinName
            End If
        Next
    End Sub

    <Extension>
    Public Sub SetSubcellularLocation(proj As GenBankProject, protein_hits As HitCollection())
        proj.annotations.transporter = protein_hits
        proj.annotations.membrane_proteins = proj.annotations.transporter _
            .Select(Function(hits) RankTerm.RankTopTerm(hits)) _
            .IteratesALL _
            .ToArray

        Dim ccIndex = proj.annotations.membrane_proteins _
            .TopTerms _
            .ToDictionary(Function(a)
                              Return a.queryName
                          End Function)

        For Each gene As GeneTable In proj.gene_table
            If ccIndex.ContainsKey(gene.locus_id) Then
                gene.subcellular_location = ccIndex(gene.locus_id).term
            End If
        Next
    End Sub

    <Extension>
    Public Sub SetTranscriptFactors(proj As GenBankProject, protein_hits As HitCollection())
        proj.annotations.tf_hits = protein_hits
        proj.annotations.transcript_factors = proj.annotations.tf_hits _
            .Select(Function(hits)
                        Return hits.AssignTFFamilyHit()
                    End Function) _
            .Where(Function(ec) Not ec Is Nothing) _
            .ToArray

        Dim tfIndex = proj.annotations.transcript_factors _
            .GroupBy(Function(a) a.QueryName) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a _
                                  .OrderByDescending(Function(i) i.identities) _
                                  .First
                          End Function)

        For Each gene As GeneTable In proj.gene_table
            If tfIndex.ContainsKey(gene.locus_id) Then
                gene.TF = tfIndex(gene.locus_id).HitName
            End If
        Next
    End Sub

    <Extension>
    Public Sub SetKEGG(proj As GenBankProject, protein_hits As HitCollection())
        proj.annotations.ko_hits = protein_hits
        proj.annotations.ko = protein_hits _
            .Select(Function(hits)
                        Return RankTerm.RankTopTerm(hits, makeTerm:=Function(a)
                                                                        Dim t = a.hitName.Split("|"c)

                                                                        If t.Length > 1 Then
                                                                            Return t(1)
                                                                        Else
                                                                            Return ""
                                                                        End If
                                                                    End Function)
                    End Function) _
            .IteratesALL _
            .ToArray

        Dim koIndex = proj.annotations.ko _
            .TopTerms _
            .ToDictionary(Function(a)
                              Return a.queryName
                          End Function)

        For Each gene As GeneTable In proj.gene_table
            If koIndex.ContainsKey(gene.locus_id) Then
                gene.KO = koIndex(gene.locus_id).term
            End If
        Next
    End Sub
End Module
