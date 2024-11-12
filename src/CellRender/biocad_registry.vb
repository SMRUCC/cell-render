Imports Oracle.LinuxCompatibility.MySQL.Uri

Public Class biocad_registry : Inherits biocad_registryModel.db_registry

    Public Sub New(mysqli As ConnectionUri)
        MyBase.New(mysqli)
    End Sub
End Class
