Imports CellBuilder
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Assembly.NCBI.GenBank
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Internal.[Object]
Imports SMRUCC.Rsharp.Runtime.Interop

<Package("project")>
<RTypeExport("genbank_project", GetType(GenBankProject))>
Module ProjectBuilder

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

End Module
