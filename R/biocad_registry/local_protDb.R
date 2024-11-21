#' Create local fasta database for protein alignment
#' 
#' @details this function will pull all protein sequence from the biocad_registry, 
#' and then build local fasta sequence database with header formats:
#' 
#' ```
#' > cad_id ec_number|function description 
#' ```
#' 
#' a placeholder symbol ``-`` will be placed after the cad_id if the ``ec_number`` 
#' of the protein is missing. for multiple ec_number, a comma seperator will be 
#' used for join the ec_number list.
#' 
const local_protDb = function(cad_registry) {

}