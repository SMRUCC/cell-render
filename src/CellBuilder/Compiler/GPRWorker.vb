Imports Microsoft.VisualBasic.ComponentModel.Collection
Imports Microsoft.VisualBasic.Language
Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Math.Scripting.MathExpression
Imports Microsoft.VisualBasic.MIME.application.json
Imports Microsoft.VisualBasic.Scripting.Expressions
Imports Microsoft.VisualBasic.Text.Xml.Models
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.genomics.GCModeller.CompilerServices
Imports SMRUCC.genomics.GCModeller.CompilerServices.GPRLink
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model
Imports SMRUCC.genomics.Interops.NCBI.Extensions.Pipeline

Public Class GPRWorker

    ReadOnly worker As MetabolicAssociator
    ReadOnly proj As GenBankProject
    ReadOnly registry As IDataRegistry

    Public Property enzyme_cutoff As Double = 450

    Sub New(proj As GenBankProject, registry As IDataRegistry)
        Dim pathways As GPRLink.Pathway() = registry.GetPathways.ToArray

        Me.proj = proj
        Me.registry = registry
        Me.worker = New MetabolicAssociator(New GPRParameters, proj.gene_table, pathways)
    End Sub

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

    Private Function FormatCompoundId(id As UInteger) As String
        Dim fullid As String = "BioCAD" & id.ToString.PadLeft(11, "0"c)
        Dim model As WebJSON.Molecule = registry.GetMoleculeDataByID(id)

        Return If(model.symbol.StringEmpty, fullid, model.symbol)
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

    Public Function CreateMetabolismNetwork(genes As Dictionary(Of String, gene)) As MetabolismStructure
        Dim scores = worker.AssociateGenesToReactions.ToArray
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
        Dim transporter As Dictionary(Of String, RankTerm) = Compiler.ProteinLocations(
            From prot As RankTerm
            In annoSet.membrane_proteins
            Where prot.term Like membranes
        )

        Call $"processing of {enzymes.Count} enzyme annotations".debug

        Dim missing_enzyme As New List(Of String)

        For Each enzyme As ECNumberAnnotation In From e As ECNumberAnnotation
                                                 In enzymes.Values
                                                 Where e.Score > enzyme_cutoff
            Dim ec_number As String = enzyme.EC
            Dim list = registry.GetAssociatedReactions(enzyme, simple:=False)

            If list Is Nothing Then
                Call missing_enzyme.Add(ec_number)
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

        If missing_enzyme.Any Then
            Call $"missing {missing_enzyme.Distinct.Count} metabolic network inside registry which is associated with enzymes: {missing_enzyme.Distinct.JoinBy(", ")}!".warning
        End If

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
            .Select(Function(id) registry.GetMoleculeDataByID(id)) _
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

End Class
