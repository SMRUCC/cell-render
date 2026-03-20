Imports System.IO
Imports System.Runtime.CompilerServices
Imports Microsoft.VisualBasic.ApplicationServices.Zip
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Serialization.JSON
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models
Imports SMRUCC.genomics.Metagenomics
Imports SMRUCC.genomics.SequenceModel.FASTA

Public Module ProjectIO

    Public Function Load(filepath As String) As GenBankProject
        Using s As Stream = filepath.Open(FileMode.Open, doClear:=False, [readOnly]:=True)
            Return Load(s)
        End Using
    End Function

    <Extension>
    Private Iterator Function LoadHitCollection(zip As ZipStream, file As String) As IEnumerable(Of HitCollection)
        Dim lines As IEnumerable(Of String) = zip.ReadLines(file)

        If lines Is Nothing Then
            Return
        Else
            Dim hit As HitCollection

            For Each line As String In lines
                If Not line.StringEmpty(, True) Then
                    hit = line.LoadJSON(Of HitCollection)(throwEx:=False)

                    If Not hit Is Nothing Then
                        Yield hit
                    End If
                End If
            Next
        End If
    End Function

    Public Function Load(s As Stream) As GenBankProject
        Using zip As New ZipStream(s, is_readonly:=True)
            Dim source_json As String = zip.ReadAllText("/source.json")
            Dim source_nt As Dictionary(Of String, Integer) = zip.ReadAllText("/source.txt").LoadJSON(Of Dictionary(Of String, Integer))
            Dim nucl_fasta As FastaSeq() = FastaFile.DocParser(zip.ReadLines("/genes.txt")).ToArray
            Dim prot_fasta As FastaSeq() = FastaFile.DocParser(zip.ReadLines("/proteins.txt")).ToArray
            Dim tss_fasta As FastaSeq() = FastaFile.DocParser(zip.ReadLines("/tss_upstream.txt")).ToArray
            Dim genes As GeneTable() = zip.ReadLines("/genes.jsonl") _
                .SafeQuery _
                .Select(Function(line) line.LoadJSON(Of GeneTable)) _
                .Where(Function(line) Not line Is Nothing) _
                .ToArray

            Dim enzyme_hits As HitCollection() = zip.LoadHitCollection("/localblast/enzyme_hits.jsonl").ToArray
            Dim operon_hits As HitCollection() = zip.LoadHitCollection("/localblast/operon_hits.jsonl").ToArray
            Dim tf_hits As HitCollection() = zip.LoadHitCollection("/localblast/tf_hits.jsonl").ToArray
            Dim transport_blast As HitCollection() = zip.LoadHitCollection("/localblast/transporter.jsonl").ToArray

            Dim operons As AnnotatedOperon() = zip.ReadLines("/localblast/operons.jsonl") _
                .SafeQuery _
                .Select(Function(line) line.LoadJSON(Of AnnotatedOperon)(throwEx:=False)) _
                .Where(Function(line) Not line Is Nothing) _
                .ToArray
            Dim ec_numbers As Dictionary(Of String, ECNumberAnnotation) = zip.ReadLines("/localblast/ec_numbers.jsonl") _
                .SafeQuery _
                .Select(Function(line) line.LoadJSON(Of ECNumberAnnotation)(throwEx:=False)) _
                .Where(Function(line) Not line Is Nothing) _
                .ToDictionary(Function(e) e.gene_id)
            Dim tfbs As MotifMatch() = zip.ReadLines("/tfbs.jsonl") _
                .SafeQuery _
                .Select(Function(line) line.LoadJSON(Of MotifMatch)(throwEx:=False)) _
                .Where(Function(line) Not line Is Nothing) _
                .ToArray
            Dim tfset As BestHit() = zip.ReadLines("/localblast/transcript_factors.jsonl") _
                .SafeQuery _
                .Select(Function(line) line.LoadJSON(Of BestHit)(throwEx:=False)) _
                .Where(Function(line) Not line Is Nothing) _
                .ToArray
            Dim membranes As RankTerm() = zip.ReadLines("/localblast/membrane_factors.jsonl") _
                .SafeQuery _
                .Select(Function(line) line.LoadJSON(Of RankTerm)(throwEx:=False)) _
                .Where(Function(line) Not line Is Nothing) _
                .ToArray

            Dim geneSeqIndex = nucl_fasta.ToDictionary(Function(a) a.Title, Function(a) a.SequenceData)
            Dim protSeqIndex = prot_fasta.ToDictionary(Function(a) a.Title, Function(a) a.SequenceData)
            Dim tssSiteIndex = tss_fasta.ToDictionary(Function(a) a.Title, Function(a) a.SequenceData)
            Dim tfbs_groups = tfbs _
                .Where(Function(a) Not a.title Is Nothing) _
                .GroupBy(Function(a) a.title) _
                .ToDictionary(Function(a) a.Key,
                                Function(a)
                                    Return a.ToArray
                                End Function)

            Return New GenBankProject With {
                .enzyme_hits = enzyme_hits,
                .nt = source_nt,
                .taxonomy = source_json.LoadJSON(Of Taxonomy)(throwEx:=False),
                .gene_table = genes,
                .genes = geneSeqIndex,
                .proteins = protSeqIndex,
                .tss_upstream = tssSiteIndex,
                .ec_numbers = ec_numbers,
                .operon_hits = operon_hits,
                .operons = operons,
                .tfbs_hits = tfbs_groups,
                .tf_hits = tf_hits,
                .transcript_factors = tfset,
                .transporter = transport_blast,
                .membrane_proteins = membranes
            }
        End Using
    End Function

    <Extension>
    Public Sub SaveZip(proj As GenBankProject, filepath As String)
        Using zip As New ZipStream(filepath.Open(FileMode.OpenOrCreate, doClear:=True, [readOnly]:=False))
            Call zip.WriteText(proj.taxonomy.GetJson, "/source.json")
            Call zip.WriteText(proj.nt.GetJson, "/source.txt")

            Call proj.DumpGeneFasta(zip.OpenFile("/genes.txt", access:=FileAccess.Write))
            Call proj.DumpProteinFasta(zip.OpenFile("/proteins.txt", access:=FileAccess.Write))
            Call proj.DumpTSSUpstreamFasta(zip.OpenFile("/tss_upstream.txt", access:=FileAccess.Write))

            Call zip.WriteLines(proj.gene_table.SafeQuery.Select(Function(a) a.GetJson), "/genes.jsonl")
            Call zip.WriteLines(proj.enzyme_hits.SafeQuery.Select(Function(q) q.GetJson), "/localblast/enzyme_hits.jsonl")
            Call zip.WriteLines(proj.operon_hits.SafeQuery.Select(Function(q) q.GetJson), "/localblast/operon_hits.jsonl")
            Call zip.WriteLines(proj.tf_hits.SafeQuery.Select(Function(q) q.GetJson), "/localblast/tf_hits.jsonl")
            Call zip.WriteLines(proj.ec_numbers.SafeQuery.Select(Function(e) e.Value.GetJson), "/localblast/ec_numbers.jsonl")
            Call zip.WriteLines(proj.transcript_factors.SafeQuery.Select(Function(e) e.GetJson), "/localblast/transcript_factors.jsonl")
            Call zip.WriteLines(proj.operons.SafeQuery.Select(Function(e) e.GetJson), "/localblast/operons.jsonl")
            Call zip.WriteLines(proj.tfbs_hits.SafeQuery.Values.IteratesALL.Select(Function(e) e.GetJson), "/tfbs.jsonl")

            Call zip.WriteLines(proj.transporter.SafeQuery.Select(Function(e) e.GetJson), "/localblast/transporter.jsonl")
            Call zip.WriteLines(proj.membrane_proteins.SafeQuery.Select(Function(e) e.GetJson), "/localblast/membrane_factors.jsonl")
        End Using
    End Sub
End Module
