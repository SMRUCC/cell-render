Imports Oracle.LinuxCompatibility.MySQL.Reflection.DbAttributes

Public Class reaction_view

    <DatabaseField> Public Property reaction_id As UInteger
    <DatabaseField> Public Property molecule_id As UInteger
    <DatabaseField> Public Property db_xref As String
    <DatabaseField> Public Property side As String
    <DatabaseField> Public Property factor As Double

End Class
