
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports SMRUCC.genomics.Analysis.HTS.DataFrame
Imports SMRUCC.Rsharp.Runtime.Interop

<Package("cad_lab")>
<RTypeExport("cad_lab", GetType(cad_lab))>
Module CadLab

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="cad_lab"></param>
    ''' <param name="exp_id">
    ''' the experiment project id
    ''' </param>
    ''' <param name="dynaimics"></param>
    ''' <returns></returns>
    <ExportAPI("save_molecules")>
    Public Function save_molecules(cad_lab As cad_lab, exp_id As String, dynaimics As Matrix)

    End Function

End Module
