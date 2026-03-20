Imports Microsoft.VisualBasic.Serialization.JSON

' Operon的类型枚举
Public Enum OperonType
    Conserved  ' 保守的
    Insertion  ' 插入突变
    Deletion   ' 缺失突变
End Enum

' 用于表示最终注释结果的Operon结构
Public Class AnnotatedOperon
    Public Property OperonID As String
    Public Property name As String
    Public Property Type As OperonType
    Public Property Genes As String() ' 组成此Operon的基因组上的基因
    Public Property Scores As Double() ' 每个基因对当前OperonID的投票总得分
    Public Property KnownGeneIds As String() ' 参考Operon中应有的基因ID
    Public Property InsertedGeneIds As String() ' 插入的新基因ID
    Public Property MissingGeneIds As String() ' 缺失的基因ID

    Public Property strand As String
    Public Property left As Integer
    Public Property right As Integer

    Public Overrides Function ToString() As String
        Return $"{Type.Description} #{OperonID} at {strand}:{left}-{right} with {Genes.TryCount} gene members {Genes.GetJson}"
    End Function
End Class

