Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Serialization.JSON
Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model

Public Class DataRepository : Implements IDataRegistry

#Region "in-memory cache data"
    ReadOnly cachedOperon As WebJSON.Operon()
    ReadOnly cachedReactions As EnzymeQuery(Of EnzymeTag)
    ReadOnly cachedMolecules As Dictionary(Of String, WebJSON.Molecule)
    ReadOnly cachedExpansion As Dictionary(Of String, WebJSON.Reaction())
#End Region

    Public ReadOnly Property HasOperonDataCache As Boolean
        Get
            Return cachedOperon IsNot Nothing
        End Get
    End Property

    Public ReadOnly Property HasReactionDataCache As Boolean
        Get
            Return cachedReactions IsNot Nothing
        End Get
    End Property

    Public ReadOnly Property HasMoleculeDataCache As Boolean
        Get
            Return cachedMolecules IsNot Nothing
        End Get
    End Property

    Public ReadOnly Property HasExpansionNetworkDataCache As Boolean
        Get
            Return cachedExpansion IsNot Nothing
        End Get
    End Property

    Dim opt As New QueryOptions

    ''' <summary>
    ''' construct a data repository without local cache data
    ''' </summary>
    Sub New()
    End Sub

    Sub New(cache_dir As String)
        Call "start to load local cahced database files".info

        Dim network = $"{cache_dir}/metabolic_network.jsonl".LoadJSONL(Of WebJSON.Reaction).ToArray

        cachedOperon = $"{cache_dir}/all_operons.json".LoadJsonFile(Of WebJSON.Operon())(throwEx:=False)
        cachedMolecules = $"{cache_dir}/molecules.jsonl".LoadJSONL(Of WebJSON.Molecule).ToDictionary(Function(m) m.id)
        cachedReactions = New EnzymeQuery(Of EnzymeTag)(TagEnzymeString(network))
        cachedExpansion = (From rxn In network Where rxn.law.IsNullOrEmpty) _
            .Select(Function(r)
                        Return r.left.JoinIterates(r.right).Select(Function(c) (c.molecule_id, r))
                    End Function) _
            .IteratesALL _
            .GroupBy(Function(c) c.molecule_id) _
            .ToDictionary(Function(c) c.Key.ToString,
                          Function(c)
                              Return c.Select(Function(i) i.r).GroupBy(Function(r) r.guid).Select(Function(r) r.First).ToArray
                          End Function)

        If cachedOperon Is Nothing Then cachedOperon = {}

        Call "load cached database from a given cache dir:".info
        Call $" * {cachedOperon.Length} known operons".info
        Call $" * {cachedReactions.Size} known enzyme reaction network".info
        Call $" * {cachedExpansion.Count} reaction network expansions".info
        Call $" * {cachedMolecules.Count} associated metabolites".info
    End Sub

    Private Shared Function TagEnzymeString(network As IEnumerable(Of WebJSON.Reaction)) As IEnumerable(Of EnzymeTag)
        Dim enzymeNetwork = From rxn In network Where Not rxn.law.IsNullOrEmpty
        Dim tagsData = enzymeNetwork.Select(Function(r) r.law.Select(Function(ec) New EnzymeTag(ec.ec_number, r))).IteratesALL

        Return tagsData
    End Function

    Public Function GetAllKnownOperons() As WebJSON.Operon()
        ' just read cache data for local test
        ' not used cache dir for save web request data
        Return cachedOperon.ToArray
    End Function

    Public Function GetMoleculeDataByID(id As UInteger) As WebJSON.Molecule Implements IDataRegistry.GetMoleculeDataByID
        Return cachedMolecules.TryGetValue(id.ToString, [default]:=cachedMolecules.TryGetValue("BioCAD" & id.ToString.PadLeft(11, "0"c)))
    End Function

    Public Function GetAssociatedReactions(enzyme As IEnzymeObject, Optional simple As Boolean = False) As Dictionary(Of String, WebJSON.Reaction) Implements IDataRegistry.GetAssociatedReactions
        Dim list As IEnumerable(Of WebJSON.Reaction) =
            From tag As EnzymeTag
            In cachedReactions.Query(enzyme.ECNumber, opt.EnzymeFuzzyMatch, opt.EnzymeMaxFuzzyLevel)
            Select tag.Reaction
            Group By Reaction.guid Into Group
            Select Group.First
        Dim asso As Dictionary(Of String, WebJSON.Reaction) = list.ToDictionary(Function(r) r.guid)

        Return asso
    End Function

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="registry_id">id of the molecule</param>
    ''' <returns></returns>
    Public Function ExpandNetworkByCompound(registry_id As String) As Dictionary(Of String, WebJSON.Reaction) Implements IDataRegistry.ExpandNetworkByCompound
        Dim list = cachedExpansion.TryGetValue(registry_id)

        If Not list Is Nothing Then
            Return list.ToDictionary(Function(r) r.guid)
        Else
            Return Nothing
        End If
    End Function

    Public Function SetOptions(opt As QueryOptions) As IDataRegistry Implements IDataRegistry.SetOptions
        Me.opt = opt
        Return Me
    End Function
End Class
