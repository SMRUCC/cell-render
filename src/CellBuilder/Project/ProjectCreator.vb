Imports Microsoft.VisualBasic.Linq
Imports SMRUCC.genomics.Assembly.NCBI.GenBank
Imports SMRUCC.genomics.Assembly.NCBI.GenBank.GBFF.Keywords.FEATURES
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.ContextModel.Promoter
Imports SMRUCC.genomics.Metagenomics
Imports SMRUCC.genomics.SequenceModel.FASTA

Public Module ProjectCreator

    Public Function FromGenBank(replicons As IEnumerable(Of GBFF.File)) As GenBankProject
        Dim nucl As New Dictionary(Of String, String)
        Dim prot As New Dictionary(Of String, String)
        Dim genes As New List(Of GeneTable)
        Dim size As New Dictionary(Of String, Integer)
        Dim tax As Taxonomy = Nothing
        Dim tss As New Dictionary(Of String, String)

        For Each replicon As GBFF.File In replicons
            Dim nt As New FastaSeq({"nt"}, replicon.Origin.SequenceData)

            If tax Is Nothing Then
                tax = replicon.Source.GetTaxonomy
            End If

            Call replicon _
                .EnumerateGeneFeatures(ORF:=False) _
                .GroupBy(Function(a) a.Query(FeatureQualifiers.locus_tag)) _
                .ToDictionary(Function(a) a.Key,
                              Function(a)
                                  Return a.First.SequenceData
                              End Function) _
                .DoCall(Sub(list) Call nucl.AddRange(list, replaceDuplicated:=True))

            Call replicon _
                .ExportProteins_Short(True) _
                .ToDictionary(Function(a) a.Title,
                              Function(a)
                                  Return a.SequenceData
                              End Function) _
                .DoCall(Sub(list) Call prot.AddRange(list, replaceDuplicated:=True))

            Dim geneSet As GeneTable() = replicon _
                .EnumerateGeneFeatures(ORF:=False) _
                .ExportTable _
                .ToArray

            Call genes.AddRange(geneSet)
            Call size.Add(replicon.Locus.AccessionID,
                          replicon.Origin.SequenceData.Length)

            Dim gene_upstreamSet = geneSet _
                .Select(Function(gene)
                            Return (gene.locus_id, gene.GetUpstreamSeq(nt, 150))
                        End Function) _
                .ToArray

            Call gene_upstreamSet _
                .GroupBy(Function(a) a.locus_id) _
                .ToDictionary(Function(a) a.Key,
                              Function(a)
                                  Return Strings.UCase(a.First.Item2.SequenceData)
                              End Function) _
                .DoCall(Sub(list)
                            Call tss.AddRange(list, replaceDuplicated:=True)
                        End Sub)
        Next

        Dim proj As New GenBankProject With {
            .taxonomy = tax,
            .nt = size,
            .genes = nucl,
            .proteins = prot,
            .gene_table = genes.ToArray,
            .tss_upstream = tss
        }

        Return proj
    End Function
End Module
