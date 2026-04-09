Imports SMRUCC.genomics.ComponentModel.Annotation
Imports SMRUCC.genomics.GCModeller.ModellingEngine.Model

Public Class EnzymeTag : Implements IEnzymeObject

    Public Property ECNumber As String Implements IEnzymeObject.ECNumber
    Public Property Reaction As WebJSON.Reaction

    Sub New(ec$, rxn As WebJSON.Reaction)
        ECNumber = ec
        Reaction = rxn
    End Sub

    Public Overrides Function ToString() As String
        Return Reaction.ToString
    End Function

End Class
