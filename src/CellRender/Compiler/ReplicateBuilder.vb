Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar
Imports Microsoft.VisualBasic.ApplicationServices.Terminal.ProgressBar.Tqdm
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular.Vector

Public Class ReplicateBuilder

    ReadOnly compiler As Compiler
    ReadOnly genes As New List(Of TranscriptUnit)
    ReadOnly rnas As New List(Of RNA)
    ReadOnly tax_id As ULong

    Public ReadOnly Property cad_registry As biocad_registry
        Get
            Return compiler.cad_registry
        End Get
    End Property

    Sub New(compiler As Compiler)
        Me.compiler = compiler
        Me.tax_id = ULong.Parse(compiler.tax_id)
    End Sub

    Private Function linkGene(gene_info As GeneTable, ByRef find As gene_molecule) As gene
        ' fetch gene information from database
        Dim findMol = cad_registry.molecule _
            .where(field("`molecule`.type") = compiler.dna_term,
                   field("xref_id") = gene_info.locus_id) _
            .find(Of biocad_registryModel.molecule)

        If findMol Is Nothing Then
            findMol = cad_registry.db_xrefs _
                .left_join("molecule") _
                .on(field("molecule.id") = field("db_xrefs.obj_id")) _
                .where(field("xref").in({
                           gene_info.locus_id,
                           gene_info.ProteinId,
                           gene_info.UniprotSwissProt,
                           gene_info.UniprotTrEMBL}, nullFilter:=True)) _
                .order_by("parent") _
                .find(Of biocad_registryModel.molecule)("`molecule`.*")
        End If

        If findMol Is Nothing Then
            Dim warn As String = $"missing gene model from the registry: {gene_info}"
            Call warn.Warning
            Call VBDebugger.EchoLine(warn)
            Return Nothing
        End If

        If findMol.parent > 0 Then
            findMol = cad_registry.molecule _
                .where(field("id") = findMol.parent) _
                .find(Of biocad_registryModel.molecule)

            If findMol Is Nothing Then
                Dim warn As String = $"missing parent replicon for the gene model: {gene_info}"
                Call warn.Warning
                Call VBDebugger.EchoLine(warn)
                Return Nothing
            End If
        End If

        find = cad_registry.molecule _
            .left_join("sequence_graph") _
            .on(field("`sequence_graph`.molecule_id") = field("`molecule`.id")) _
            .where(field("`molecule`.id") = findMol.id) _
            .find(Of gene_molecule)("`molecule`.id", "xref_id", "name", "note", "sequence")

        ' missing current gene item inside database
        If find Is Nothing Then
            Dim warn As String = $"missing gene model from the registry: {gene_info}"

            Call warn.Warning
            Call VBDebugger.EchoLine(warn)

            Return Nothing
        End If

        Dim rna = RNAComposition _
            .FromNtSequence(find.sequence, gene_info.locus_id) _
            .CreateVector
        Dim find_prot = cad_registry.molecule _
            .left_join("sequence_graph") _
            .on(field("`sequence_graph`.molecule_id") = field("`molecule`.id")) _
            .where(field("`molecule`.parent") = find.id) _
            .find(Of gene_molecule)("`molecule`.id", "molecule.xref_id", "sequence")
        Dim gene As New gene(gene_info.Location) With {
            .locus_tag = gene_info.locus_id,
            .product = find.note,
            .nucleotide_base = rna
        }

        If find_prot Is Nothing AndAlso Not gene_info.ProteinId.StringEmpty(, True) Then
            find_prot = cad_registry.db_xrefs _
                .left_join("sequence_graph") _
                .on(field("`sequence_graph`.molecule_id") = field("obj_id")) _
                .where(field("type") = compiler.polypeptide_term, field("xref") = gene_info.ProteinId) _
                .find(Of gene_molecule)("molecule_id AS id", "xref AS xref_id", "sequence")
        End If

        If Not find_prot Is Nothing Then
            ' find a protein sequnece
            ' is CDS/ORF
            gene.protein_id = find_prot.id
            gene.amino_acid = ProteinComposition _
                .FromRefSeq(find_prot.sequence, find_prot.xref_id) _
                .CreateVector
            gene.type = RNATypes.mRNA
        Else
            gene.type = RNAType(gene_info.type)
            ' no protein sequence could be found
            ' is rRNA or tRNA or other kind of RNA
            Call rnas.Add(New RNA With {
                .gene = gene.locus_tag,
                .type = gene.type,
                .val = gene_info.commonName
            })
        End If

        Return gene
    End Function

    ''' <summary>
    ''' contains CDS/tRNA/rRNA
    ''' </summary>
    ''' <returns></returns>
    Public Function BuildGenome() As replicon
        Dim bar As Tqdm.ProgressBar = Nothing
        Dim operons = cad_registry.conserved_cluster _
            .where(field("tax_id") = tax_id) _
            .select(Of biocad_registryModel.conserved_cluster)
        Dim tu_units As New Index(Of String)
        Dim template_index = compiler.template.ToDictionary(Function(a) a.locus_id)

        Call VBDebugger.EchoLine("compile of the genome model, pull gene and proteins.")
        Call VBDebugger.EchoLine("processing of the conserved transcript unit")

        For Each operon As biocad_registryModel.conserved_cluster In TqdmWrapper.Wrap(operons, bar:=bar)
            Dim tu = cad_registry.cluster_link _
                .left_join("molecule") _
                .on(field("`molecule`.id") = field("gene_id")) _
                .where(field("cluster_id") = operon.id) _
                .select(Of biocad_registryModel.molecule)("molecule.*")
            ' matches from the templates list
            Dim members As New List(Of gene)

            Call bar.SetLabel(operon.db_xref)

            For Each unit_gene In tu
                Dim gene_info As GeneTable = template_index.TryGetValue(unit_gene.xref_id)

                If Not gene_info Is Nothing Then
                    Dim find As gene_molecule = Nothing
                    Dim gene_model As gene = linkGene(gene_info, find)

                    If gene_model Is Nothing OrElse find Is Nothing Then
                        Continue For
                    End If

                    Call members.Add(gene_model)
                    Call template_index.Remove(unit_gene.xref_id)
                End If
            Next

            If tu.Length <> members.Count Then
                Call $"missing gene object for operon {operon.name}".Warning
            End If

            If members.Count > 0 Then
                Call genes.Add(New TranscriptUnit With {
                      .id = operon.db_xref,
                      .name = operon.name,
                      .genes = members.ToArray
                })
            End If
        Next

        Call VBDebugger.EchoLine($"create {genes.Count} conserved operons!")
        Call VBDebugger.EchoLine("processing other genes that not inside the conserved operons")

        For Each gene_info As GeneTable In TqdmWrapper.Wrap(template_index.Values, bar:=bar)
            Call bar.SetLabel(gene_info.ToString)

            Dim find As gene_molecule = Nothing
            Dim gene_model As gene = linkGene(gene_info, find)

            If gene_model Is Nothing OrElse find Is Nothing Then
                Continue For
            End If

            Dim gene_tu As New TranscriptUnit With {
                .id = find.id,
                .name = gene_info.geneName,
                .genes = {gene_model}
            }

            Call genes.Add(gene_tu)
        Next

        Return New replicon With {
            .genomeName = "",
            .isPlasmid = False,
            .operons = genes.ToArray,
            .RNAs = rnas.ToArray
        }
    End Function

    Private Shared Function RNAType(s As String) As RNATypes
        Select Case Strings.Trim(s).ToLower
            ' broken data!
            Case "cds" : Return RNATypes.micsRNA
            Case "trna"
                Return RNATypes.tRNA
            Case "rrna"
                Return RNATypes.ribosomalRNA
            Case Else
                Return RNATypes.micsRNA
        End Select
    End Function
End Class
