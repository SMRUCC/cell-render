Imports System.Runtime.CompilerServices
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports Oracle.LinuxCompatibility.MySQL.Uri

Public Class biocad_registry : Inherits biocad_registryModel.db_registry

    Public ReadOnly Property molecule As Model
        Get
            Return m_molecule
        End Get
    End Property

    Public ReadOnly Property regulation_graph As Model
        Get
            Return m_regulation_graph
        End Get
    End Property

    Public Sub New(mysqli As ConnectionUri)
        MyBase.New(mysqli)
    End Sub

    <MethodImpl(MethodImplOptions.AggressiveInlining)>
    Public Function GetVocabulary(term As String) As biocad_registryModel.vocabulary
        Return m_vocabulary.where(field("term") = term).find(Of biocad_registryModel.vocabulary)
    End Function
End Class
