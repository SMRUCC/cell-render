Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.WebJSON

Public Interface IDataRegistry

    Function GetMoleculeDataByID(id As UInteger) As Molecule
    Function GetAssociatedReactions(ec_number As String, Optional simple As Boolean = False) As Dictionary(Of String, Reaction)

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="registry_id">id of the molecule</param>
    ''' <returns></returns>
    Function ExpandNetworkByCompound(registry_id As String) As Dictionary(Of String, Reaction)

End Interface
