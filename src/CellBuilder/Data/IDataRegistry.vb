Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.CompilerServices.GPRLink
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
    Function GetPathways() As IEnumerable(Of Pathway)

End Interface

Public Class QueryOptions

    Public Property EnzymeFuzzyMatch As Boolean = False
    ''' <summary>
    ''' max allow fuzzy level, value should be c(1,2,3,4), which means the class levels of thee ec_number(which has 4 ranks).
    ''' 
    ''' example as 4, means only allows fuzzy matches on the last digit
    ''' 3, means allows fuzzy matches on the level 3 and level 4 digits
    ''' </summary>
    ''' <returns></returns>
    Public Property EnzymeMaxFuzzyLevel As Integer = 4

End Class