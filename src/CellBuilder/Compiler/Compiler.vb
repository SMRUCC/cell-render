Imports Microsoft.VisualBasic.CommandLine
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.MIME.application.json
Imports Microsoft.VisualBasic.Text.Xml.Models
Imports SMRUCC.genomics.Analysis.metaTraits.Traitar
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline

Public Class Compiler : Inherits Compiler(Of VirtualCell)

    Friend ReadOnly proj As GenBankProject
    Friend ReadOnly registry As IDataRegistry
    Friend ReadOnly motifSites As Dictionary(Of String, MotifMatch())
    Friend ReadOnly defaultName As String

    Public Property enzyme_cutoff As Double = 450

    Sub New(proj As GenBankProject, registry As IDataRegistry, Optional defaultName As String = Nothing)
        Dim annoSet As AnnotationSet = proj.annotations

        Me.defaultName = defaultName
        Me.registry = registry
        Me.motifSites = annoSet.tfbs_hits.Values _
            .IteratesALL _
            .Where(Function(a) a.pvalue < 0.05) _
            .GroupBy(Function(a) ParsePWMID(a)) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a.ToArray
                          End Function)
        Me.proj = proj
    End Sub

    Private Shared Function ParsePWMID(tfbs As MotifMatch) As String
        Dim id As String() = tfbs.seeds(0).Split
        Dim pwm_id As String = id(1)
        Return pwm_id
    End Function

    Protected Overrides Function PreCompile(args As CommandLine) As Integer
        Dim name As String = args("--name") Or defaultName

        If name.StringEmpty Then
            If proj.taxonomy Is Nothing Then
                name = Now.ToString.MD5
            Else
                name = proj.taxonomy.scientificName.Replace(" "c, "_").StringReplace("[-\._]{2,}", "_")
            End If
        End If

        m_compiledModel = New VirtualCell With {
            .taxonomy = proj.taxonomy,
            .properties = New SMRUCC.genomics.GCModeller.CompilerServices.[Property],
            .cellular_id = name
        }

        Call $"target genome taxonomy information: {proj.taxonomy.GetJson}".info
        Call $"cell name: {m_compiledModel.cellular_id}".info

        Return 0
    End Function

    Protected Overrides Function CompileImpl(args As CommandLine) As Integer
        Dim genomics As New TRNWorker(Me)

        m_compiledModel.genome = genomics.BuildGenome()
        m_compiledModel.metabolismStructure = CreateMetabolismNetwork(m_compiledModel.genome)
        m_compiledModel.traits = New Traits With {
            .phenotype = proj.annotations.traits _
                .Where(Function(p) p.status = "trained" AndAlso p.confidence >= 0.85) _
                .Select(Function(p)
                            Return p.ToPhenotype
                        End Function) _
                .ToArray
        }

        Call "link the cellular component success!".info
        Call "compile virtual cell model job done!".info

        Return 0
    End Function

    Private Function CreateMetabolismNetwork(genome As Genome) As MetabolismStructure
        Dim geneSet As Dictionary(Of String, gene) = genome.replicons _
            .Select(Function(r) r.GetGeneList) _
            .IteratesALL _
            .GroupBy(Function(a) a.locus_tag) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a.First
                          End Function)
        Dim gpr As New GPRWorker(proj, registry) With {
            .enzyme_cutoff = enzyme_cutoff
        }

        Return gpr.CreateMetabolismNetwork(geneSet)
    End Function

    Friend Shared Function ProteinLocations(list As IEnumerable(Of RankTerm)) As Dictionary(Of String, RankTerm)
        Return list _
            .GroupBy(Function(a) a.queryName) _
            .Select(Function(a)
                        Return a.OrderByDescending(Function(i) i.score).First
                    End Function) _
            .ToDictionary(Function(t)
                              Return t.queryName
                          End Function)
    End Function
End Class
