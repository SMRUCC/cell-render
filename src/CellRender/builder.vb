Imports System.IO
Imports CellBuilder
Imports Microsoft.VisualBasic.CommandLine
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Assembly.NCBI.GenBank
Imports SMRUCC.genomics.GCModeller.Assembly.GCMarkupLanguage.v2
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Components
Imports SMRUCC.Rsharp.Runtime.Internal.[Object]
Imports SMRUCC.Rsharp.Runtime.Interop
Imports RInternal = SMRUCC.Rsharp.Runtime.Internal

<Package("project")>
<RTypeExport("genbank_project", GetType(GenBankProject))>
Module ProjectBuilder

    Sub Main()
        Call RInternal.generic.add("writeBin", GetType(GenBankProject), AddressOf saveProject)
    End Sub

    Private Function saveProject(proj As GenBankProject, args As list, env As Environment) As Object
        Return save(proj, args!con, env)
    End Function

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="replicons">a vector of the ncbi genbank object of the genome replicons</param>
    ''' <param name="env"></param>
    ''' <returns></returns>
    <ExportAPI("new")>
    <RApiReturn(GetType(GenBankProject))>
    Public Function create(<RRawVectorArgument> replicons As Object, Optional env As Environment = Nothing) As Object
        Dim pull As pipeline = pipeline.TryCreatePipeline(Of GBFF.File)(replicons, env)

        If pull.isError Then
            Return pull.getError
        End If

        Return New ProjectCreator().FromGenBank(pull.populates(Of GBFF.File)(env))
    End Function

    <ExportAPI("save")>
    Public Function save(proj As GenBankProject, file As Object, Optional env As Environment = Nothing) As Object
        Dim is_file As Boolean = False
        Dim s = SMRUCC.Rsharp.GetFileStream(file, IO.FileAccess.Write, env, is_filepath:=is_file)

        If s Like GetType(Message) Then
            Return s.TryCast(Of Message)
        Else
            Call ProjectIO.SaveZip(proj, s.TryCast(Of Stream))

            Try
                Call s.TryCast(Of Stream).Flush()
            Catch ex As Exception

            End Try
        End If

        Try
            If is_file Then
                Call s.TryCast(Of Stream).Dispose()
            End If
        Catch ex As Exception

        End Try

        Return Nothing
    End Function

    <ExportAPI("load")>
    Public Function load(<RRawVectorArgument> file As Object, Optional env As Environment = Nothing) As Object
        Dim is_file As Boolean = False
        Dim s = SMRUCC.Rsharp.GetFileStream(file, IO.FileAccess.Read, env, is_filepath:=is_file)

        If s Like GetType(Message) Then
            Return s.TryCast(Of Message)
        Else
            Return ProjectIO.Load(s.TryCast(Of Stream))
        End If
    End Function

    <ExportAPI("scientific_name")>
    Public Function taxonomy_name(proj As GenBankProject) As String
        Return proj.taxonomy.scientificName
    End Function

    ''' <summary>
    ''' build the genbank project as a virtual cell model 
    ''' </summary>
    ''' <param name="proj"></param>
    ''' <param name="datapool"></param>
    ''' <param name="vcell_name"></param>
    ''' <returns></returns>
    <ExportAPI("build")>
    Public Function build(proj As GenBankProject, datapool As DataRepository, Optional vcell_name As String = Nothing) As VirtualCell
        Dim compiler As New Compiler(proj, datapool, defaultName:=vcell_name)
        Dim args As CommandLine = CommandLine.Parse($"--build --name ""{vcell_name}""")
        Dim model As VirtualCell = compiler.Compile(args)

        Return model
    End Function
End Module
