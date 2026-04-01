Imports Microsoft.VisualBasic.Linq
Imports Microsoft.VisualBasic.Serialization.JSON
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model

Public Class DataRepository : Implements IDataRegistry

#Region "in-memory cache data"
    ReadOnly cachedOperon As WebJSON.Operon()
    ReadOnly cachedReactions As Dictionary(Of String, WebJSON.Reaction())
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
        cachedReactions = (From rxn In network Where Not rxn.law.IsNullOrEmpty) _
            .Select(Function(r) r.law.Select(Function(ec) (ec.ec_number, r))) _
            .IteratesALL _
            .GroupBy(Function(r) r.ec_number) _
            .ToDictionary(Function(r) r.Key,
                          Function(r)
                              Return r.Select(Function(i) i.r).ToArray
                          End Function)
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
        Call $" * {cachedReactions.Count} known enzyme reaction network".info
        Call $" * {cachedExpansion.Count} reaction network expansions".info
        Call $" * {cachedMolecules.Count} associated metabolites".info
    End Sub

    Public Function GetAllKnownOperons() As WebJSON.Operon()
        ' just read cache data for local test
        ' not used cache dir for save web request data
        Return cachedOperon.ToArray
    End Function

    Public Function GetMoleculeDataByID(id As UInteger) As WebJSON.Molecule Implements IDataRegistry.GetMoleculeDataByID
        Return cachedMolecules.TryGetValue(id.ToString, [default]:=cachedMolecules.TryGetValue("BioCAD" & id.ToString.PadLeft(11, "0"c)))
    End Function

    Public Function GetAssociatedReactions(ec_number As String, Optional simple As Boolean = False) As Dictionary(Of String, WebJSON.Reaction) Implements IDataRegistry.GetAssociatedReactions
        Dim list = cachedReactions.TryGetValue(ec_number)

        If Not list Is Nothing Then
            Return list.ToDictionary(Function(r) r.guid)
        Else
            Return Nothing
        End If
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
End Class
