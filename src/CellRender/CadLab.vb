
Imports Microsoft.VisualBasic.CommandLine.Reflection
Imports Microsoft.VisualBasic.Scripting.MetaData
Imports Microsoft.VisualBasic.Serialization.BinaryDumping
Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports SMRUCC.genomics.Analysis.HTS.DataFrame
Imports SMRUCC.Rsharp.Runtime
Imports SMRUCC.Rsharp.Runtime.Interop
Imports RInternal = SMRUCC.Rsharp.Runtime.Internal

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
    Public Function save_molecules(cad_lab As cad_lab, exp_id As String, dynaimics As Matrix, Optional env As Environment = Nothing) As Object
        Dim exp = cad_lab.experiment _
            .where(field("proj_id") = exp_id) _
            .find(Of biocad_labModel.experiment)

        If exp Is Nothing Then
            cad_lab.experiment.add(
                field("proj_id") = exp_id,
                field("name") = exp_id
            )

            exp = cad_lab.experiment _
                .where(field("proj_id") = exp_id) _
                .order_by("id", desc:=True) _
                .find(Of biocad_labModel.experiment)

            If exp Is Nothing Then
                Return RInternal.debug.stop($"create experiment project with id '{exp_id}' error!", env)
            End If
        End If

        Dim tagdata = dynaimics.tag.GetTagValue(":")
        Dim cell = cad_lab.virtualcell _
            .where(field("run_id") = exp.id,
                   field("name") = tagdata.Value,
                   field("taxid") = tagdata.Name) _
            .find(Of biocad_labModel.virtualcell)

        If cell Is Nothing Then
            cad_lab.virtualcell.add(
                field("run_id") = exp.id,
                field("name") = tagdata.Value,
                field("taxid") = tagdata.Name
            )

            cell = cad_lab.virtualcell _
                .where(field("run_id") = exp.id,
                       field("name") = tagdata.Value,
                       field("taxid") = tagdata.Name) _
                .order_by("id", desc:=True) _
                .find(Of biocad_labModel.virtualcell)

            If cell Is Nothing Then
                Return RInternal.debug.stop($"create virtual cell with name '{tagdata.Value}' and taxid '{tagdata.Name}' error!", env)
            End If
        End If

        Static encoder As New NetworkByteOrderBuffer

        For Each mol As DataFrameRow In dynaimics.expression
            Dim moldata = cad_lab.dynamics _
                .where(field("cell_id") = cell.id,
                       field("mol_id") = mol.geneID) _
                .find(Of biocad_labModel.dynamics)

            If moldata Is Nothing Then
                cad_lab.dynamics.add(
                    field("cell_id") = cell.id,
                    field("mol_id") = mol.geneID,
                    field("x0") = mol(Scan0),
                    field("dynamics") = encoder.Base64String(mol.experiments, gzip:=True)
                )
            End If
        Next

        Return Nothing
    End Function

End Module
