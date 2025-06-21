Imports CellRender.biocad_labModel
Imports Oracle.LinuxCompatibility.MySQL
Imports Oracle.LinuxCompatibility.MySQL.Uri

Public Class cad_lab : Inherits biocad_labModel.db_cad_lab

    Public ReadOnly Property dynamics As TableModel(Of dynamics)
        Get
            Return m_dynamics
        End Get
    End Property

    Public ReadOnly Property experiment As TableModel(Of experiment)
        Get
            Return m_experiment
        End Get
    End Property

    Public ReadOnly Property virtualcell As TableModel(Of virtualcell)
        Get
            Return m_virtualcell
        End Get
    End Property

    Public Sub New(mysqli As ConnectionUri)
        MyBase.New(mysqli)
    End Sub
End Class
