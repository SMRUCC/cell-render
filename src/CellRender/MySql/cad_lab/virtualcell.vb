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
''' DROP TABLE IF EXISTS `virtualcell`;
''' /*!40101 SET @saved_cs_client     = @@character_set_client */;
''' /*!50503 SET character_set_client = utf8mb4 */;
''' CREATE TABLE `virtualcell` (
'''   `id` int unsigned NOT NULL AUTO_INCREMENT,
'''   `run_id` int unsigned NOT NULL COMMENT 'the experiment run id',
'''   `name` varchar(255) NOT NULL,
'''   `taxid` int NOT NULL DEFAULT '0',
'''   `add_time` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
'''   PRIMARY KEY (`id`),
'''   UNIQUE KEY `id_UNIQUE` (`id`),
'''   KEY `project_info_idx` (`run_id`)
''' ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;
''' /*!40101 SET character_set_client = @saved_cs_client */;
''' 
''' --
''' -- Dumping events for database 'cad_lab'
''' --
''' 
''' --
''' -- Dumping routines for database 'cad_lab'
''' --
''' /*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;
''' 
''' /*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
''' /*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
''' /*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
''' /*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
''' /*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
''' /*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
''' /*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
''' 
''' -- Dump completed on 2025-06-21 15:36:11
''' ```
''' </summary>
''' <remarks></remarks>
<Oracle.LinuxCompatibility.MySQL.Reflection.DbAttributes.TableName("virtualcell", Database:="cad_lab", SchemaSQL:="
CREATE TABLE `virtualcell` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `run_id` int unsigned NOT NULL COMMENT 'the experiment run id',
  `name` varchar(255) NOT NULL,
  `taxid` int NOT NULL DEFAULT '0',
  `add_time` datetime NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  UNIQUE KEY `id_UNIQUE` (`id`),
  KEY `project_info_idx` (`run_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;")>
Public Class virtualcell: Inherits Oracle.LinuxCompatibility.MySQL.MySQLTable
#Region "Public Property Mapping To Database Fields"
    <DatabaseField("id"), PrimaryKey, AutoIncrement, NotNull, DataType(MySqlDbType.UInt32, "11"), Column(Name:="id"), XmlAttribute> Public Property id As UInteger
''' <summary>
''' the experiment run id
''' </summary>
''' <value></value>
''' <returns></returns>
''' <remarks></remarks>
    <DatabaseField("run_id"), NotNull, DataType(MySqlDbType.UInt32, "11"), Column(Name:="run_id")> Public Property run_id As UInteger
    <DatabaseField("name"), NotNull, DataType(MySqlDbType.VarChar, "255"), Column(Name:="name")> Public Property name As String
    <DatabaseField("taxid"), NotNull, DataType(MySqlDbType.Int32, "11"), Column(Name:="taxid")> Public Property taxid As Integer
    <DatabaseField("add_time"), NotNull, DataType(MySqlDbType.DateTime), Column(Name:="add_time")> Public Property add_time As Date
#End Region
#Region "Public SQL Interface"
#Region "Interface SQL"
    Friend Shared ReadOnly INSERT_SQL$ = 
        <SQL>INSERT INTO `virtualcell` (`run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}');</SQL>

    Friend Shared ReadOnly INSERT_AI_SQL$ = 
        <SQL>INSERT INTO `virtualcell` (`id`, `run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');</SQL>

    Friend Shared ReadOnly REPLACE_SQL$ = 
        <SQL>REPLACE INTO `virtualcell` (`run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}');</SQL>

    Friend Shared ReadOnly REPLACE_AI_SQL$ = 
        <SQL>REPLACE INTO `virtualcell` (`id`, `run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');</SQL>

    Friend Shared ReadOnly DELETE_SQL$ =
        <SQL>DELETE FROM `virtualcell` WHERE `id` = '{0}';</SQL>

    Friend Shared ReadOnly UPDATE_SQL$ = 
        <SQL>UPDATE `virtualcell` SET `id`='{0}', `run_id`='{1}', `name`='{2}', `taxid`='{3}', `add_time`='{4}' WHERE `id` = '{5}';</SQL>

#End Region

''' <summary>
''' ```SQL
''' DELETE FROM `virtualcell` WHERE `id` = '{0}';
''' ```
''' </summary>
    Public Overrides Function GetDeleteSQL() As String
        Return String.Format(DELETE_SQL, id)
    End Function

''' <summary>
''' ```SQL
''' INSERT INTO `virtualcell` (`id`, `run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetInsertSQL() As String
        Return String.Format(INSERT_SQL, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time))
    End Function

''' <summary>
''' ```SQL
''' INSERT INTO `virtualcell` (`id`, `run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetInsertSQL(AI As Boolean) As String
        If AI Then
        Return String.Format(INSERT_AI_SQL, id, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time))
        Else
        Return String.Format(INSERT_SQL, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time))
        End If
    End Function

''' <summary>
''' <see cref="GetInsertSQL"/>
''' </summary>
    Public Overrides Function GetDumpInsertValue(AI As Boolean) As String
        If AI Then
            Return $"('{id}', '{run_id}', '{name}', '{taxid}', '{add_time.ToString("yyyy-MM-dd hh:mm:ss")}')"
        Else
            Return $"('{run_id}', '{name}', '{taxid}', '{add_time.ToString("yyyy-MM-dd hh:mm:ss")}')"
        End If
    End Function


''' <summary>
''' ```SQL
''' REPLACE INTO `virtualcell` (`id`, `run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetReplaceSQL() As String
        Return String.Format(REPLACE_SQL, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time))
    End Function

''' <summary>
''' ```SQL
''' REPLACE INTO `virtualcell` (`id`, `run_id`, `name`, `taxid`, `add_time`) VALUES ('{0}', '{1}', '{2}', '{3}', '{4}');
''' ```
''' </summary>
    Public Overrides Function GetReplaceSQL(AI As Boolean) As String
        If AI Then
        Return String.Format(REPLACE_AI_SQL, id, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time))
        Else
        Return String.Format(REPLACE_SQL, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time))
        End If
    End Function

''' <summary>
''' ```SQL
''' UPDATE `virtualcell` SET `id`='{0}', `run_id`='{1}', `name`='{2}', `taxid`='{3}', `add_time`='{4}' WHERE `id` = '{5}';
''' ```
''' </summary>
    Public Overrides Function GetUpdateSQL() As String
        Return String.Format(UPDATE_SQL, id, run_id, name, taxid, MySqlScript.ToMySqlDateTimeString(add_time), id)
    End Function
#End Region

''' <summary>
                     ''' Memberwise clone of current table Object.
                     ''' </summary>
                     Public Function Clone() As virtualcell
                         Return DirectCast(MyClass.MemberwiseClone, virtualcell)
                     End Function
End Class


End Namespace
