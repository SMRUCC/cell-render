Imports CADRegistry
Imports CellBuilder
Imports Microsoft.VisualBasic.CommandLine
Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Language
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Math.Scripting.MathExpression
Imports Microsoft.VisualBasic.MIME.application.json
Imports Microsoft.VisualBasic.Text.Xml.Models
Imports SMRUCC.genomics.Analysis.SequenceTools.SequencePatterns
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.Cellular.Vector
Imports SMRUCC.genomics.Interops.NCBI.Extensions.LocalBLAST.Application.BBH
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline
Imports SMRUCC.genomics.Model.OperonMapper
Imports SMRUCC.genomics.SequenceModel.NucleotideModels.Translation

Public Class Compiler : Inherits Compiler(Of VirtualCell)

    ReadOnly proj As GenBankProject
    ReadOnly registry As RegistryUrl
    ReadOnly motifSites As Dictionary(Of String, MotifMatch())
    ReadOnly defaultName As String

    Public Property enzyme_cutoff As Double = 450

    Sub New(proj As GenBankProject, Optional serverUrl As String = RegistryUrl.defaultServer, Optional defaultName As String = Nothing)
        Dim annoSet As AnnotationSet = proj.annotations

        Me.defaultName = defaultName

        If serverUrl.ToLower.StartsWith("http://") OrElse serverUrl.ToLower.StartsWith("https://") Then
            Me.registry = New RegistryUrl(serverUrl)
        Else
            Me.registry = New RegistryUrl(RegistryUrl.defaultServer, serverUrl)
        End If

        Me.motifSites = annoSet.tfbs_hits.Values _
            .IteratesALL _
            .Where(Function(a) a.identities > 0.97) _
            .GroupBy(Function(a) a.seeds(0)) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a.ToArray
                          End Function)
        Me.proj = proj
    End Sub

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
        m_compiledModel.genome = BuildGenome()
        m_compiledModel.metabolismStructure = m_compiledModel.genome.replicons _
            .Select(Function(r) r.GetGeneList) _
            .IteratesALL _
            .GroupBy(Function(a) a.locus_tag) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a.First
                          End Function) _
            .DoCall(Function(list)
                        Return CreateMetabolismNetwork(list)
                    End Function)

        Call "link the cellular component success!".info
        Call "compile virtual cell model job done!".info

        Return 0
    End Function

    Private Iterator Function BuildLaws(reaction As WebJSON.Reaction, enzyme As ECNumberAnnotation, modelProteinId As String) As IEnumerable(Of Catalysis)
        For Each law As WebJSON.LawData In reaction.law.SafeQuery
            Dim pars = law.params.Keys.ToArray
            Dim args As KineticsParameter() = law.params _
                .Select(Function(a)
                            Return CreateParameter(a, modelProteinId)
                        End Function) _
                .ToArray

            Yield New Catalysis With {
                .reaction = reaction.guid,
                .temperature = 36,
                .PH = 7.0,
                .formula = New FunctionElement With {
                    .lambda = law.lambda,
                    .name = enzyme.EC,
                    .parameters = pars
                },
                .parameter = args
            }
        Next
    End Function

    Private Function CreateParameter(a As KeyValuePair(Of String, String), modelProteinId As String) As KineticsParameter
        If a.Value.IsNumeric Then
            Return New KineticsParameter With {
                .name = a.Key,
                .value = Val(a.Value),
                .isModifier = False
            }
        ElseIf a.Value.StartsWith("ENZ_") Then
            Return New KineticsParameter With {
                .name = a.Key,
                .value = 0,
                .isModifier = False,
                .target = modelProteinId
            }
        ElseIf a.Value.IsPattern("BioCAD\d+") Then
            Dim m As String = FormatCompoundId(UInteger.Parse(a.Value.Match("\d+")))
            Dim k As New KineticsParameter With {
                .name = a.Key,
                .value = 0,
                .isModifier = False,
                .target = m
            }

            Return k
        Else
            Return New KineticsParameter With {
                .name = a.Key,
                .value = 0,
                .isModifier = False,
                .target = a.Value
            }
        End If
    End Function

    Private Function CreateMetabolismNetwork(genes As Dictionary(Of String, gene)) As MetabolismStructure
        Dim annoSet As AnnotationSet = proj.annotations
        Dim enzymes As Dictionary(Of String, ECNumberAnnotation) = annoSet.ec_numbers
        Dim network As New Dictionary(Of String, WebJSON.Reaction)
        Dim ec_numbers As New Dictionary(Of String, List(Of String))
        Dim enzymeModels As New List(Of Enzyme)
        Dim geneIndex = proj.gene_table _
            .GroupBy(Function(a) a.locus_id) _
            .ToDictionary(Function(a) a.Key,
                          Function(a)
                              Return a.First
                          End Function)

        Static membranes As Index(Of String) = {"Cell_inner_membrane", "Cell_membrane", "Cell_outer_membrane"}

        Dim membraneTransport As New List(Of (ECNumberAnnotation, String, String()))
        Dim transporter As Dictionary(Of String, RankTerm) = ProteinLocations(
            From prot As RankTerm
            In annoSet.membrane_proteins
            Where prot.term Like membranes
        )

        Call $"processing of {enzymes.Count} enzyme annotations".debug

        For Each enzyme As ECNumberAnnotation In From e As ECNumberAnnotation
                                                 In enzymes.Values
                                                 Where e.Score > enzyme_cutoff
            Dim ec_number As String = enzyme.EC
            Dim list = registry.GetAssociatedReactions(enzyme.EC, simple:=False)

            If list Is Nothing Then
                Call $"missing metabolic network inside registry which is associated with enzyme {enzyme.EC}!".warning
                Continue For
            Else
                Dim gene As GeneTable = geneIndex(enzyme.gene_id)
                Dim translate_id As String = If(gene.ProteinId, gene.locus_id & "_translate")
                Dim modelProteinId As String = "Protein[" & translate_id & "]"
                Dim model As New Enzyme With {
                    .ECNumber = enzyme.EC,
                    .proteinID = modelProteinId,
                    .catalysis = list.Values _
                         .Select(Function(reaction) BuildLaws(reaction, enzyme, modelProteinId)) _
                         .IteratesALL _
                         .GroupBy(Function(a) a.GetJson.MD5) _
                         .Select(Function(a) a.First) _
                         .ToArray
                }

                If transporter.ContainsKey(gene.locus_id) Then
                    Call membraneTransport.Add((enzyme, transporter(gene.locus_id).term, list.Keys.ToArray))
                End If

                Call enzymeModels.Add(model)
            End If

            Call network.AddRange(From r
                                  In list
                                  Where r.Value.left _
                                      .JoinIterates(r.Value.right) _
                                      .All(Function(a) a.molecule_id > 0), replaceDuplicated:=True)

            For Each guid As String In list.Keys
                If Not ec_numbers.ContainsKey(guid) Then
                    Call ec_numbers.Add(guid, New List(Of String))
                End If

                ec_numbers(guid).Add(ec_number)
            Next
        Next

        Call $"load {network.Count} enzymatic reactions!".debug

        Dim none_enzymatic = ExpandNetwork(network).ToArray
        Dim metabolites As Compound() = CreateCompoundModel(network, none_enzymatic) _
            .OrderBy(Function(c) c.ID) _
            .ToArray

        Return New MetabolismStructure With {
            .compounds = metabolites,
            .reactions = New ReactionGroup With {
                .enzymatic = CreateEnzymaticNetwork(network, ec_numbers).ToArray,
                .transportation = membraneTransport _
                    .Select(Function(a)
                                Return a.Item3.Select(Function(rxn_id) (rxn_id, a.Item1.gene_id, cc:=a.Item2))
                            End Function) _
                    .IteratesALL _
                    .GroupBy(Function(a) a.rxn_id) _
                    .Select(Function(a)
                                Return New Transportation With {
                                    .guid = a.Key,
                                    .enzymes = a.Select(Function(i) i.gene_id).Distinct.ToArray,
                                    .membrane = a.Select(Function(i) i.cc).Distinct.ToArray
                                }
                            End Function) _
                    .ToArray,
                .none_enzymatic = none_enzymatic
            },
            .enzymes = enzymeModels.ToArray
        }
    End Function

    Private Iterator Function ExpandNetwork(network As Dictionary(Of String, WebJSON.Reaction)) As IEnumerable(Of Reaction)
        Dim compounds_id As UInteger() = network.Values _
           .Select(Function(r) r.left.JoinIterates(r.right)) _
           .IteratesALL _
           .GroupBy(Function(a) a.molecule_id) _
           .Keys
        Dim pending As New Queue(Of UInteger)(compounds_id)
        Dim cache As New Dictionary(Of UInteger, WebJSON.Reaction())

        Call "start to expends the reaction network...".debug

        Do While pending.Count > 0
            Dim mol_id As UInteger = pending.Dequeue

            If cache.ContainsKey(mol_id) Then
                Continue Do
            End If

            Dim expansions As Dictionary(Of String, WebJSON.Reaction) = If(
                registry.ExpandNetworkByCompound(mol_id),
                New Dictionary(Of String, WebJSON.Reaction)
            )

            Call cache.Add(mol_id, expansions.Values.ToArray)

            For Each r As WebJSON.Reaction In expansions.Values
                Dim new_compounds As UInteger() = r.left.JoinIterates(r.right) _
                    .Select(Function(c) c.molecule_id) _
                    .ToArray

                For Each id As UInteger In new_compounds
                    Call pending.Enqueue(id)
                Next
            Next
        Loop

        For Each groupdata In cache.Values.IteratesALL.GroupBy(Function(a) a.guid)
            Dim reaction As WebJSON.Reaction = groupdata.First
            Dim hash_id As String = reaction.guid
            Dim left = MakeSubstrates(reaction.left).ToArray
            Dim right = MakeSubstrates(reaction.right).ToArray
            Dim model As New Reaction With {
                .bounds = {5, 5},
                .compartment = Nothing,
                .ec_number = Nothing,
                .ID = hash_id,
                .is_enzymatic = False,
                .name = reaction.name,
                .note = reaction.reaction,
                .substrate = left,
                .product = right
            }

            Yield model
        Next
    End Function

    Private Iterator Function MakeSubstrates(list As IEnumerable(Of WebJSON.Substrate)) As IEnumerable(Of CompoundFactor)
        For Each c As WebJSON.Substrate In list
            Yield New CompoundFactor With {
                .factor = c.factor,
                .compound = FormatCompoundId(c.molecule_id),
                .cid = c.molecule_id
            }
        Next
    End Function

    Private Iterator Function CreateEnzymaticNetwork(network As Dictionary(Of String, WebJSON.Reaction), ec_numbers As Dictionary(Of String, List(Of String))) As IEnumerable(Of Reaction)
        For Each a As WebJSON.Reaction In network.Values
            Dim hash_id As String = a.guid
            Dim left = MakeSubstrates(a.left).ToArray
            Dim right = MakeSubstrates(a.right).ToArray

            Yield New Reaction With {
                .ID = hash_id,
                .ec_number = ec_numbers(a.guid).Distinct.ToArray,
                .bounds = {5, 5},
                .is_enzymatic = True,
                .name = a.name,
                .note = a.reaction,
                .substrate = left,
                .product = right
            }
        Next

        Call "create the enzymatic network success".info
    End Function

    Private Function ExtractCompoundModelId(network As Dictionary(Of String, WebJSON.Reaction), none_enzymatic As Reaction()) As IEnumerable(Of UInteger)
        Dim cset1 As UInteger() = network.Values _
            .Select(Function(r) r.left.JoinIterates(r.right)) _
            .IteratesALL _
            .GroupBy(Function(a) a.molecule_id) _
            .Keys
        Dim cset2 As UInteger() = none_enzymatic _
            .Select(Function(r)
                        Return r.substrate.JoinIterates(r.product)
                    End Function) _
            .IteratesALL _
            .Select(Function(f) f.cid) _
            .ToArray

        Return cset1.JoinIterates(cset2).Distinct
    End Function

    Private Iterator Function CreateCompoundModel(network As Dictionary(Of String, WebJSON.Reaction), none_enzymatic As Reaction()) As IEnumerable(Of Compound)
        Dim compounds_id As UInteger() = ExtractCompoundModelId(network, none_enzymatic).ToArray
        Dim metadata As WebJSON.Molecule() = compounds_id _
            .Select(Function(id) registry.GetMoleculeDataById(id)) _
            .Where(Function(c) Not c Is Nothing) _
            .ToArray
        Dim refs As Index(Of String) = {"BioCyc", "MetaCyc", "KEGG"}

        Call $"found {compounds_id.Length} associated metabolites!".debug

        For Each c As WebJSON.Molecule In metadata
            Dim biocyc_id As WebJSON.DBXref() = c.db_xrefs _
                .SafeQuery _
                .Where(Function(r) r.dbname Like refs) _
                .ToArray

            Yield New v2.Compound With {
                .db_xrefs = c.db_xrefs _
                    .SafeQuery _
                    .Select(Function(a) a.xref_id) _
                    .Distinct _
                    .ToArray,
                .ID = FormatCompoundId(UInteger.Parse(c.id.Match("\d+"))),
                .name = c.name,
                .referenceIds = biocyc_id _
                    .Select(Function(xi) xi.xref_id) _
                    .ToArray,
                .formula = c.formula
            }
        Next
    End Function

    Private Shared Function ProteinLocations(list As IEnumerable(Of RankTerm)) As Dictionary(Of String, RankTerm)
        Return list _
            .GroupBy(Function(a) a.queryName) _
            .Select(Function(a)
                        Return a.OrderByDescending(Function(i) i.score).First
                    End Function) _
            .ToDictionary(Function(t)
                              Return t.queryName
                          End Function)
    End Function

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="rnas">其他的RNA列表</param>
    ''' <param name="proteins"></param>
    ''' <param name="regulations"></param>
    ''' <returns></returns>
    Private Iterator Function GeneObjects(rnas As List(Of (rid$, RNA)), proteins As List(Of protein), regulations As List(Of transcription)) As IEnumerable(Of (rid$, gene))
        Dim nt As Dictionary(Of String, String) = proj.genes
        Dim RNA As RNA
        Dim annoSet As AnnotationSet = proj.annotations
        Dim tfs As Dictionary(Of String, BestHit()) = annoSet.transcript_factors _
            .GroupBy(Function(tf) tf.QueryName) _
            .ToDictionary(Function(t) t.Key,
                          Function(t)
                              Return t.ToArray
                          End Function)
        Dim protein_id As String
        Dim transporter As Dictionary(Of String, RankTerm) = ProteinLocations(annoSet.membrane_proteins)

        Call $"processing compile of {nt.Count} genes!".debug

        For Each gene As GeneTable In proj.gene_table
            Dim nt_seq As String = nt(gene.locus_id)
            Dim bases As NumericVector = RNAComposition.FromNtSequence(nt_seq, gene.locus_id & "_rna").CreateVector
            Dim residues As NumericVector = Nothing
            Dim gene_type As RNATypes
            Dim translate_id As String = If(gene.ProteinId, gene.locus_id & "_translate")
            Dim isTF As Boolean = tfs.ContainsKey(gene.locus_id)

            If Not gene.translation.StringEmpty Then
                residues = ProteinComposition.FromRefSeq(gene.translation, translate_id).CreateVector
                gene_type = RNATypes.mRNA
                protein_id = "Protein[" & translate_id & "]"

                If isTF Then
                    Call regulations.AddRange(RegulationNetwork(protein_id, gene.locus_id, annotation:=tfs(gene.locus_id)))
                End If

                Call proteins.Add(New protein With {
                    .name = protein_id,
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
                            .name = protein_id,
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
    End Function

    Private Function BuildGenome() As Genome
        Dim RNAs As New List(Of (rid$, RNA))
        Dim proteins As New List(Of protein)
        Dim regulationNetwork As New List(Of transcription)
        Dim replicons As New List(Of replicon)
        Dim annoSet As AnnotationSet = proj.annotations

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
            If motifSites.ContainsKey(hit.HitName) Then
                For Each site As MotifMatch In motifSites(hit.HitName)
                    Yield New transcription With {
                        .regulator = regulator,
                        .targets = {site.title},
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
            End If
        Next
    End Function

    Private Function FormatCompoundId(id As UInteger) As String
        Dim fullid As String = "BioCAD" & id.ToString.PadLeft(11, "0"c)
        Dim model As WebJSON.Molecule = registry.GetMoleculeDataById(id)

        Return If(model.symbol.StringEmpty, fullid, model.symbol)
    End Function
End Class
