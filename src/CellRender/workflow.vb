Imports System.IO
Imports CellBuilder
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.Assembly.KEGG.DBGET.bGetObject
Imports SMRUCC.genomics.Assembly.KEGG.WebServices.XML
Imports SMRUCC.genomics.Assembly.MetaCyc.Schema
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Tasks.Models
Imports SMRUCC.genomics.SequenceModel.FASTA
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Components
Imports SMRUCC.Rsharp.Runtime.Internal.[Object]
Imports SMRUCC.Rsharp.Runtime.Interop
Imports RInternal = SMRUCC.Rsharp.Runtime.Internal

''' <summary>
''' annotation workflow
''' </summary>
<Package("workflow")>
Module workflow

    ''' <summary>
    ''' extract of the tss upstream location site sequence data
    ''' </summary>
    ''' <param name="proj"></param>
    ''' <returns></returns>
    <ExportAPI("tss_upstream")>
    <RApiReturn(GetType(FastaSeq))>
    Public Function tss_upstream(proj As GenBankProject) As Object
        Return proj _
            .DumpTSSUpstreamFasta _
            .ToArray
    End Function

    ''' <summary>
    ''' extract of the protein fasta sequence data to file
    ''' </summary>
    ''' <param name="proj"></param>
    ''' <param name="file"></param>
    ''' <param name="env"></param>
    ''' <returns></returns>
    <ExportAPI("save_proteins")>
    Public Function save_proteins(proj As GenBankProject, file As Object, Optional env As Environment = Nothing) As Object
        Dim s = SMRUCC.Rsharp.GetFileStream(file, FileAccess.Write, env)

        If s Like GetType(Message) Then
            Return s.TryCast(Of Message)
        End If

        Call proj.DumpProteinFasta(s.TryCast(Of Stream))

        Return Nothing
    End Function

    ''' <summary>
    ''' get enzyme annotation result table from the project model
    ''' </summary>
    ''' <param name="proj"></param>
    ''' <returns></returns>
    <ExportAPI("enzyme_table")>
    Public Function enzyme_table(proj As GenBankProject) As ECNumberAnnotation()
        If proj Is Nothing OrElse proj.annotations Is Nothing Then
            Return Nothing
        Else
            Return proj.annotations.ec_numbers.Values.ToArray
        End If
    End Function

    <ExportAPI("set_blastp_result")>
    Public Function set_blastp_result(<RRawVectorArgument> blastp_hits As Object, proj As GenBankProject, group As String, Optional env As Environment = Nothing) As Object
        Dim pull As pipeline = pipeline.TryCreatePipeline(Of HitCollection)(blastp_hits, env)

        If pull.isError Then
            Return pull.getError
        ElseIf proj.annotations Is Nothing Then
            proj.annotations = New AnnotationSet
        End If

        Dim protein_hits As HitCollection() = pull.populates(Of HitCollection)(env).ToArray

        Select Case Strings.LCase(group)
            Case "ec_number"
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
                    End If
                Next

            Case "subcellular_location"
                proj.annotations.transporter = protein_hits
                proj.annotations.membrane_proteins = proj.annotations.transporter _
                    .Select(Function(hits) RankTerm.RankTopTerm(hits)) _
                    .IteratesALL _
                    .ToArray

            Case "transcript_factor"
                proj.annotations.tf_hits = protein_hits
                proj.annotations.transcript_factors = proj.annotations.tf_hits _
                    .Select(Function(hits)
                                Return hits.AssignTFFamilyHit()
                            End Function) _
                    .Where(Function(ec) Not ec Is Nothing) _
                    .ToArray

            Case Else
                Return RInternal.debug.stop($"unknown annotation group of '{group}'", env)
        End Select

        Return Nothing
    End Function

    <ExportAPI("set_tfbs")>
    Public Function set_tfbs(proj As GenBankProject, <RRawVectorArgument> tfbs As Object, Optional env As Environment = Nothing) As Object
        Dim pull As pipeline = pipeline.TryCreatePipeline(Of MotifMatch)(tfbs, env)

        If pull.isError Then
            Return pull.getError
        ElseIf proj.annotations Is Nothing Then
            proj.annotations = New AnnotationSet
        End If

        proj.annotations.tfbs_hits = pull _
            .populates(Of MotifMatch)(env) _
            .GroupBy(Function(m) m.title) _
            .ToDictionary(Function(m) m.Key,
                          Function(m)
                              Return m.ToArray
                          End Function)

        Return proj
    End Function

    <ExportAPI("open_datapool")>
    Public Function open_datapool(dir As String, Optional enzyme_fuzzy As Boolean = False) As DataRepository
        Return New DataRepository(dir).SetOptions(New QueryOptions With {.EnzymeFuzzyMatch = enzyme_fuzzy, .EnzymeMaxFuzzyLevel = 4})
    End Function

    <ExportAPI("set_kegg_pathways")>
    Public Function set_kegg_pathways(repo As DataRepository,
                                      <RRawVectorArgument(GetType(Map))> maps As Object,
                                      <RRawVectorArgument(GetType(Reaction))> reactions As Object,
                                      Optional env As Environment = Nothing) As Object

        Dim pullMaps As pipeline = pipeline.TryCreatePipeline(Of Map)(maps, env)
        Dim pullReactions As pipeline = pipeline.TryCreatePipeline(Of Reaction)(reactions, env)

        If pullMaps.isError Then
            Return pullMaps.getError
        ElseIf pullReactions.isError Then
            Return pullReactions.getError
        End If

        Static filterMaps As Index(Of String) = {
            "map01100",  ' Metabolic pathways
            "map01110",  ' Biosynthesis of secondary metabolites
            "map01120"   ' Microbial metabolism in diverse environments
        }

        Dim cleanMaps As IEnumerable(Of Map) = From map As Map
                                               In pullMaps.populates(Of Map)(env)
                                               Where Not map.EntryId Like filterMaps

        Return repo.SetPathwayData(GPRLink.Pathway.FromKEGGPathways(
            pathways:=cleanMaps,
            reactions:=pullReactions.populates(Of Reaction)(env))
        )
    End Function
End Module
