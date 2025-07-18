REM  Oracle.LinuxCompatibility.MySQL.CodeSolution.VisualBasic.CodeGenerator
REM  MYSQL Schema Mapper
REM      for Microsoft VisualBasic.NET 1.0.0.0

REM  Dump @6/21/2025 3:39:04 PM


Imports System.Xml.Serialization
Imports Microsoft.VisualBasic.ComponentModel.DataSourceModel.SchemaMaps
Imports Oracle.LinuxCompatibility.MySQL.Reflection.DbAttributes
Imports MySqlScript = Oracle.LinuxCompatibility.MySQL.Scripting.Extensions

Namespace biocad_labModel

''' <summary>
''' ```SQL
''' 
''' --
''' 
''' DROP TABLE IF EXISTS `dynamics`;
''' /*!40101 SET @saved_cs_client     = @@character_set_client */;
''' /*!50503 SET character_set_client = utf8mb4 */;
''' CREATE TABLE `dynamics` (
'''   `id` int unsigned NOT NULL AUTO_INCREMENT,
'''   `cell_id` int unsigned NOT NULL,
'''   `mol_id` int unsigned NOT NULL COMMENT 'the biocad registry molecule id',
'''   `dynamics` longtext NOT NULL COMMENT 'gzip compressed based 64 encoded molecule content double[] vector',
'''   `add_time` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
'''   PRIMARY KEY (`id`),
'''   UNIQUE KEY `id_UNIQUE` (`id`),
'''   KEY `cell_source_idx` (`cell_id`)
''' ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
''' /*!40101 SET character_set_client = @saved_cs_client */;
''' 
''' --
''' ```
''' </summary>
''' <remarks></remarks>
<Oracle.LinuxCompatibility.MySQL.Reflection.DbAttributes.TableName("dynamics", Database:="cad_lab", SchemaSQL:="
CREATE TABLE `dynamics` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `cell_id` int unsigned NOT NULL,
  `mol_id` int unsigned NOT NULL COMMENT 'the biocad registry molecule id',
  `dynamics` longtext NOT NULL COMMENT 'gzip compressed based 64 encoded molecule content double[] vector',
  `add_time` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `id_UNIQUE` (`id`),
  KEY `cell_source_idx` (`cell_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;")>
Public Class dynamics: Inherits Oracle.LinuxCompatibility.MySQL.MySQLTable
#Region "Public Property Mapping To Database Fields"
    <DatabaseField("id"), PrimaryKey, AutoIncrement, NotNull, DataType(MySqlDbType.UInt32, "11"), Column(Name:="id"), XmlAttribute> Public Property id As UInteger
    <DatabaseField("cell_id"), NotNull, DataType(MySqlDbType.UInt32, "11"), Column(Name:="cell_id")> Public Property cell_id As UInteger
''' <summary>
''' the biocad registry molecule id
''' </summary>
''' <value></value>
''' <returns></returns>
''' <remarks></remarks>
    <DatabaseField("mol_id"), NotNull, DataType(MySqlDbType.UInt32, "11"), Column(Name:="mol_id")> Public Property mol_id As UInteger
''' <summary>
''' gzip compressed based 64 encoded molecule content double[] vector
''' </summary>
''' <value></value>
''' <returns></returns>
''' <remarks></remarks>
    <DatabaseField("dynamics"), NotNull, DataType(MySqlDbType.Text), Column(Name:="dynamics")> Public Property dynamics As String
    <DatabaseField("add_time"), NotNull, DataType(MySqlDbType.DateTime), Column(Name:="add_time")> Public Property add_time As Date
#End Region
#Region "Public SQL Interface"
#Region "Interface SQL"
    Friend Shared ReadOnly INSERT_SQL$ = 
        <SQL>INSERT INTO `dynamics` (`cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}');</SQL>

    Friend Shared ReadOnly INSERT_AI_SQL$ = 
        <SQL>INSERT INTO `dynamics` (`id`, `cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');</SQL>

    Friend Shared ReadOnly REPLACE_SQL$ = 
        <SQL>REPLACE INTO `dynamics` (`cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}');</SQL>

    Friend Shared ReadOnly REPLACE_AI_SQL$ = 
        <SQL>REPLACE INTO `dynamics` (`id`, `cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');</SQL>

    Friend Shared ReadOnly DELETE_SQL$ =
        <SQL>DELETE FROM `dynamics` WHERE `id` = '{0}';</SQL>

    Friend Shared ReadOnly UPDATE_SQL$ = 
        <SQL>UPDATE `dynamics` SET `id`='{0}', `cell_id`='{1}', `mol_id`='{2}', `dynamics`='{3}', `add_time`='{4}' WHERE `id` = '{5}';</SQL>

#End Region

''' <summary>
''' ```SQL
''' DELETE FROM `dynamics` WHERE `id` = '{0}';
''' ```
''' </summary>
    Public Overrides Function GetDeleteSQL() As String
        Return String.Format(DELETE_SQL, id)
    End Function

''' <summary>
''' ```SQL
''' INSERT INTO `dynamics` (`id`, `cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetInsertSQL() As String
        Return String.Format(INSERT_SQL, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time))
    End Function

''' <summary>
''' ```SQL
''' INSERT INTO `dynamics` (`id`, `cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetInsertSQL(AI As Boolean) As String
        If AI Then
        Return String.Format(INSERT_AI_SQL, id, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time))
        Else
        Return String.Format(INSERT_SQL, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time))
        End If
    End Function

''' <summary>
''' <see cref="GetInsertSQL"/>
''' </summary>
    Public Overrides Function GetDumpInsertValue(AI As Boolean) As String
        If AI Then
            Return $"('{id}', '{cell_id}', '{mol_id}', '{dynamics}', '{add_time.ToString("yyyy-MM-dd hh:mm:ss")}')"
        Else
            Return $"('{cell_id}', '{mol_id}', '{dynamics}', '{add_time.ToString("yyyy-MM-dd hh:mm:ss")}')"
        End If
    End Function


''' <summary>
''' ```SQL
''' REPLACE INTO `dynamics` (`id`, `cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetReplaceSQL() As String
        Return String.Format(REPLACE_SQL, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time))
    End Function

''' <summary>
''' ```SQL
''' REPLACE INTO `dynamics` (`id`, `cell_id`, `mol_id`, `dynamics`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetReplaceSQL(AI As Boolean) As String
        If AI Then
        Return String.Format(REPLACE_AI_SQL, id, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time))
        Else
        Return String.Format(REPLACE_SQL, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time))
        End If
    End Function

''' <summary>
''' ```SQL
''' UPDATE `dynamics` SET `id`='{0}', `cell_id`='{1}', `mol_id`='{2}', `dynamics`='{3}', `add_time`='{4}' WHERE `id` = '{5}';
''' ```
''' </summary>
    Public Overrides Function GetUpdateSQL() As String
        Return String.Format(UPDATE_SQL, id, cell_id, mol_id, dynamics, MySqlScript.ToMySqlDateTimeString(add_time), id)
    End Function
#End Region

''' <summary>
                     ''' Memberwise clone of current table Object.
                     ''' </summary>
                     Public Function Clone() As dynamics
                         Return DirectCast(MyClass.MemberwiseClone, dynamics)
                     End Function
End Class


End Namespace
