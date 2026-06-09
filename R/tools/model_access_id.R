#' Extract the Model Genomic Accession ID from Replicons
#'
#' Retrieves the primary genomic accession ID from a list of GenBank replicons.
#' It prioritizes the chromosomal sequence (non-plasmid) over plasmid sequences.
#'
#' @details
#' The function iterates through the provided \code{replicons}. As soon as it encounters
#' a replicon that is not flagged as a plasmid (typically representing the main
#' chromosome), it immediately returns its accession ID.
#' If the GenBank file contains \emph{only} plasmid data, the function falls back to
#' returning the accession ID of the last plasmid processed.
#' If the input list is empty, it returns \code{NULL}.
#'
#' @param replicons A list or vector of GenBank replicon objects (e.g., parsed
#'   from a multi-replicon GenBank file).
#'
#' @return A character string representing the accession ID, or \code{NULL} if no
#'   replicons are provided.
#'
#' @seealso \code{\link[GenBank]{accession_id}} for extracting accession IDs,
#'   and \code{\link{is.plasmid}} for checking replicon type,
#'   \code{\link{make_genbank_proj_file}} which uses this function to
#'   determine model directory names in batch mode.
#'
#' @examples
#' \dontrun{
#' # Assuming 'gb_file_replicons' is a list of parsed GenBank replicon objects
#' model_id <- model_accession_id(gb_file_replicons)
#' print(model_id)
#' }
#'
#' @export
const model_accession_id = function(replicons) {
    let accession_id = NULL;

    for(rep in replicons) {
        if (!is.plasmid(rep)) {
            return(GenBank::accession_id(rep));
        } 

        accession_id <- GenBank::accession_id(rep);
    }

    # all data inside this genbank file is plasmid???
    return(accession_id);
}