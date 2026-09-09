Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Language
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Scripting.Expressions
Imports Microsoft.VisualBasic.Text.Xml.Models
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular.Vector
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Model.OperonMapper
Imports SMRUCC.genomics.SequenceModel.NucleotideModels.Translation

Public Class TRNWorker

    ReadOnly compiler As Compiler

    Sub New(compiler As Compiler)
        Me.compiler = compiler
    End Sub

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="rnas">其他的RNA列表</param>
    ''' <param name="proteins"></param>
    ''' <param name="regulations"></param>
    ''' <returns></returns>
    Private Iterator Function GeneObjects(rnas As List(Of (rid$, RNA)), proteins As List(Of protein), regulations As List(Of transcription)) As IEnumerable(Of (rid$, gene))
        Dim nt As Dictionary(Of String, String) = compiler.proj.genes
        Dim RNA As RNA
        Dim annoSet As AnnotationSet = compiler.proj.annotations
        Dim tfs As Dictionary(Of String, BestHit()) = annoSet.transcript_factors _
            .GroupBy(Function(tf) tf.QueryName) _
            .ToDictionary(Function(t) t.Key,
                          Function(t)
                              Return t.ToArray
                          End Function)
        Dim protein_id As String
        Dim transporter As Dictionary(Of String, RankTerm) = Compiler.ProteinLocations(annoSet.membrane_proteins)
        Dim missing_seqs As New List(Of String)
        Dim ec_numbers = annoSet.ec_numbers

        Call $"processing compile of {nt.Count} genes!".debug

        For Each gene As GeneTable In compiler.proj.gene_table
            Dim nt_seq As String = If(nt.TryGetValue(gene.locus_id), "")
            Dim bases As NumericVector = RNAComposition.FromNtSequence(nt_seq, gene.locus_id & "_rna").CreateVector
            Dim residues As NumericVector = Nothing
            Dim gene_type As RNATypes
            Dim translate_id As String = If(gene.ProteinId, gene.locus_id & "_translate")
            Dim isTF As Boolean = tfs.ContainsKey(gene.locus_id)
            Dim ec As ECNumberAnnotation = ec_numbers.TryGetValue(gene.locus_id)

            If nt_seq = "" Then
                Call missing_seqs.Add(nt_seq)
            End If

            If Not gene.translation.StringEmpty Then
                residues = ProteinComposition.FromRefSeq(gene.translation, translate_id).CreateVector
                gene_type = RNATypes.mRNA
                protein_id = "Protein[" & translate_id & "]"

                If isTF Then
                    Call regulations.AddRange(RegulationNetwork(protein_id, gene.locus_id, annotation:=tfs(gene.locus_id)))
                End If

                Call proteins.Add(New protein With {
                    .name = If(ec?.proteinName, protein_id),
                    .peptide_chains = {translate_id},
                    .protein_id = protein_id,
                    .cellular_location = If(transporter.ContainsKey(gene.locus_id), transporter(gene.locus_id).term, "Cytoplasm")
                })
            Else
                Select Case gene.type
                    Case "CDS"
                        Dim trans As String = TranslationTable.Translate(nt_seq)

                        residues = ProteinComposition.FromRefSeq(trans, translate_id).CreateVector
                        gene_type = RNATypes.mRNA
                        protein_id = "Protein[" & translate_id & "]"

                        If isTF Then
                            Call regulations.AddRange(RegulationNetwork(protein_id, gene.locus_id, annotation:=tfs(gene.locus_id)))
                        End If

                        Call proteins.Add(New protein With {
                            .name = If(ec?.proteinName, protein_id),
                            .peptide_chains = {translate_id},
                            .protein_id = protein_id,
                            .cellular_location = If(transporter.ContainsKey(gene.locus_id), transporter(gene.locus_id).term, "Cytoplasm")
                        })
                    Case "rRNA"
                        Dim rRNA = Strings.Trim(gene.commonName) _
                            .Replace("ribosomal RNA", "") _
                            .Trim _
                            .Split(" "c, "-"c) _
                            .FirstOrDefault

                        rRNA = Strings.Trim(rRNA).ToLower
                        gene_type = If(rRNA = "", RNATypes.micsRNA, RNATypes.ribosomalRNA)
                        RNA = New RNA With {
                            .gene = gene.locus_id,
                            .id = If(gene.commonName, gene.locus_id & "-micsRNA"),
                            .note = gene.commonName,
                            .type = gene_type,
                            .val = rRNA
                        }
                        rnas.Add((gene.replicon_accessionID, RNA))
                    Case "tRNA"
                        gene_type = If(gene.commonName.StringEmpty, RNATypes.micsRNA, RNATypes.tRNA)
                        RNA = New RNA With {
                            .gene = gene.locus_id,
                            .id = If(gene.commonName, gene.locus_id & "-micsRNA"),
                            .note = gene.commonName,
                            .type = gene_type,
                            .val = Strings.Trim(gene.commonName).Split("-"c).Last
                        }
                        rnas.Add((gene.replicon_accessionID, RNA))
                    Case Else
                        gene_type = RNATypes.micsRNA
                        RNA = New RNA With {
                            .gene = gene.locus_id,
                            .id = If(gene.commonName, gene.locus_id & "_micsRNA"),
                            .note = gene.commonName,
                            .type = gene_type,
                            .val = ""
                        }
                        rnas.Add((gene.replicon_accessionID, RNA))

                        If isTF Then
                            Call regulations.AddRange(RegulationNetwork(RNA.id, gene.locus_id, annotation:=tfs(gene.locus_id)))
                        End If
                End Select
            End If

            Dim model As New gene With {
                .locus_tag = gene.locus_id,
                .left = gene.left,
                .right = gene.right,
                .strand = gene.strand,
                .product = gene.commonName,
                .type = gene_type,
                .amino_acid = residues,
                .nucleotide_base = bases,
                .protein_id = If(residues Is Nothing, Nothing, {residues.name})
            }

            Yield (gene.replicon_accessionID, model)
        Next

        Call $"found {rnas.Count} RNA models!".debug

        If missing_seqs.Any Then
            Call $"missing {missing_seqs.Count} nucleotide sequences(missing={missing_seqs.JoinBy(", ")})!".warning
        End If
    End Function

    Public Function BuildGenome() As Genome
        Dim RNAs As New List(Of (rid$, RNA))
        Dim proteins As New List(Of protein)
        Dim regulationNetwork As New List(Of transcription)
        Dim replicons As New List(Of replicon)
        Dim annoSet As AnnotationSet = compiler.proj.annotations

        For Each replicon_group As (rid As String, genes As gene()) In GeneObjects(RNAs, proteins, regulationNetwork) _
            .GroupBy(Function(a) a.rid) _
            .Select(Function(g)
                        Return (rid:=g.Key, genes:=g.Select(Function(a) a.Item2).ToArray)
                    End Function)

            Dim geneSet As Dictionary(Of String, gene) = replicon_group.genes _
                .GroupBy(Function(a) a.locus_tag) _
                .ToDictionary(Function(a) a.Key,
                              Function(a)
                                  Return a.First
                              End Function)
            Dim operons As New List(Of TranscriptUnit)

            For Each op As AnnotatedOperon In annoSet.operons.SafeQuery
                Dim geneList = op.Genes _
                    .Where(Function(gene_id) geneSet.ContainsKey(gene_id)) _
                    .Select(Function(gene_id)
                                Return geneSet(gene_id)
                            End Function) _
                    .ToArray

                If geneList.Any Then
                    Call operons.Add(New TranscriptUnit With {
                        .id = op.OperonID,
                        .name = op.name,
                        .note = op.Type.Description,
                        .genes = geneList.ToArray
                    })
                End If
            Next

            Call $"get {operons.Count} operons was annotated in genome replicon model(id={replicon_group.rid})!".info

            Dim operon_genes As Index(Of String) = operons _
                .Select(Function(op) op.genes) _
                .IteratesALL _
                .Select(Function(gene) gene.locus_tag) _
                .Indexing

            For Each gene As gene In geneSet.Values
                If Not gene.locus_tag Like operon_genes Then
                    Call operons.Add(New TranscriptUnit(gene))
                End If
            Next

            Dim genomics As New replicon With {
                .genomeName = replicon_group.rid,
                .isPlasmid = False,
                .operons = operons.ToArray,
                .RNAs = RNAs.Where(Function(a) a.rid = .genomeName) _
                    .Select(Function(r) r.Item2) _
                    .ToArray
            }

            Call replicons.Add(genomics)
        Next

        Return New Genome With {
            .replicons = replicons.ToArray,
            .proteins = proteins.ToArray,
            .regulations = regulationNetwork.ToArray
        }
    End Function

    Private Iterator Function RegulationNetwork(regulator$, gene_id$, annotation As BestHit()) As IEnumerable(Of transcription)
        For Each hit As BestHit In annotation
            Dim pwm_id As String = hit.HitName.Split("|"c).Last

            If Not compiler.motifSites.ContainsKey(pwm_id) Then
                Continue For
            End If

            For Each site As MotifMatch In compiler.motifSites(pwm_id)
                Yield New transcription With {
                    .regulator = regulator,
                    .operonId = site.title,
                    .motif = New Motif With {
                        .left = site.start,
                        .right = site.ends,
                        .strand = "?",
                        .sequence = site.segment
                    },
                    .mode = "+"c,
                    .note = $"TFBS match: {hit.HitName}, motif score: [{site.score1:F2}, {site.score2:F2}], identities: {site.identities:P2}",
                    .centralDogma = {site.title}
                }
            Next
        Next
    End Function
End Class
