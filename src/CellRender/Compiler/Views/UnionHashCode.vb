Imports Oracle.LinuxCompatibility.MySQL.MySqlBuilder
Imports Oracle.LinuxCompatibility.MySQL.Reflection.DbAttributes

Public Class UnionHashCode

    ''' <summary>
    ''' the reacttion has the same ec_number, substrate, products if there hashcode identical
    ''' </summary>
    ''' <returns></returns>
    <DatabaseField> Public Property hashcode As String
    <DatabaseField> Public Property reaction_id As String

    Public Shared Function LoadUniqueHashCodes(cad_registry As biocad_registry) As UnionHashCode()
        Dim reaction_term As UInteger = cad_registry.getVocabulary("Reaction", "Entity Type")
        Dim hash_query = cad_registry.hashcode _
            .where(field("hashcode") <> "",
                   field("type_id") = reaction_term) _
            .group_by("hashcode") _
            .select(Of UnionHashCode)("hashcode", "GROUP_CONCAT(DISTINCT obj_id) AS reaction_id")

        Return hash_query
    End Function

End Class
