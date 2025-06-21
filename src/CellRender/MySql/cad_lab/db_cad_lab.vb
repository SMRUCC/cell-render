Imports Oracle.LinuxCompatibility.MySQL
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports Oracle.LinuxCompatibility.MySQL.Uri

Namespace biocad_labModel

Public MustInherit Class db_cad_lab : Inherits IDatabase
Protected ReadOnly m_dynamics As TableModel(Of dynamics)
Protected ReadOnly m_experiment As TableModel(Of experiment)
Protected ReadOnly m_virtualcell As TableModel(Of virtualcell)
Protected Sub New(mysqli As ConnectionUri)
Call MyBase.New(mysqli)

Me.m_dynamics = model(Of dynamics)()
Me.m_experiment = model(Of experiment)()
Me.m_virtualcell = model(Of virtualcell)()
End Sub
End Class

End Namespace
