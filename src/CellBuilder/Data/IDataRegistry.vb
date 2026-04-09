Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model.WebJSON

Public Interface IDataRegistry

    Function SetOptions(opt As QueryOptions) As IDataRegistry

    Function GetMoleculeDataByID(id As UInteger) As Molecule
    Function GetAssociatedReactions(enzyme As IEnzymeObject, Optional simple As Boolean = False) As Dictionary(Of String, Reaction)

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="registry_id">id of the molecule</param>
    ''' <returns></returns>
    Function ExpandNetworkByCompound(registry_id As String) As Dictionary(Of String, Reaction)

End Interface

Public Class QueryOptions

    Public Property EnzymeFuzzyMatch As Boolean = False

End Class