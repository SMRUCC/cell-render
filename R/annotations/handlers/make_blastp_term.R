#' Parse DIAMOND BLASTP Results into Project Session
#'
#' Reads standard DIAMOND BLASTP tabular output (m8 format) from a specified
#' working directory, groups the hits, and writes them as annotation terms into
#' a virtual cell modeling project session file. These terms are essential for
#' constructing distinct cellular networks in later steps.
#'
#' Specifically, it parses three files and maps them to the following networks:
#' \itemize{
#'   \item \code{ec_number.m8}: Extracts EC numbers for building the
#'     \strong{metabolic network}.
#'   \item \code{subcellular.m8}: Extracts subcellular localization terms for
#'     building the \strong{transmembrane network}.
#'   \item \code{transcript_factor.m8}: Extracts transcription factor hits for
#'     building the \strong{gene expression transcription regulation network}.
#' }
#'
#' @param proj_file Character. The file path to the GenBank annotation project
#'   file (\code{.gcproj}) into which the parsed terms will be written.
#' @param model_dir Character. The directory path containing the DIAMOND BLASTP
#'   output files (\code{ec_number.m8}, \code{subcellular.m8},
#'   \code{transcript_factor.m8}).
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item Injects grouped annotation terms into the project session.
#'     \item Exports an enzyme annotation table as \code{enzymes.csv} in the
#'       same directory as \code{proj_file}.
#'     \item Saves the updated project file to disk.
#'   }
#'
#' @details
#' The function loads the project file, reads each m8 file, groups the hits
#' via \code{diamond_hitgroups}, and stores them using \code{set_blastp_result}.
#' After all three annotation types are processed, it generates an enzyme
#' summary table and persists the updated project.
#'
#' @seealso \code{\link{make_terms}} which calls this function,
#'   \code{\link{make_diamond_hits}} which produces the input m8 files.
#'
#' @examples
#' \dontrun{
#' # Parse BLASTP results for a single model
#' make_blastp_term(
#'   proj_file = "release/GCF_000123/builder.gcproj",
#'   model_dir = "workspace/make_diamond_hits/GCF_000123"
#' )
#' }
#'
#' @keywords internal
#' @export
const make_blastp_term = function(proj_file, model_dir) {
    # read the diamond blastp result files from the model worker dir
    let ec_number = read_m8(file.path(model_dir, "ec_number.m8"));
    let subcellular = read_m8(file.path(model_dir, "subcellular.m8"));
    let tf_list = read_m8(file.path(model_dir, "transcript_factor.m8"));
    let proj = project::load(proj_file);

    let write_proj_session = function(blastp_result, key) {
        # write the blastp annotation terms result
        # into the current project session
        blastp_result
        |> diamond_hitgroups 
        |> set_blastp_result(proj, key)
        ;
    }

    # ec number terms for build the metabolic network
    write_proj_session(ec_number, "ec_number");
    # subcellular location terms for build the transmembrane network
    write_proj_session(subcellular, "subcellular_location");
    # tf list term hits for build the gene expression transcription regulation network
    write_proj_session(tf_list, "transcript_factor");

    let enzymes = enzyme_table(proj);
    let enzyme_file = file.path(dirname(proj_file), "enzymes.csv");

    # export enzyme annotation table file
    write.csv(enzymes, file = enzyme_file, silent = TRUE);
    # save project file
    project::save(proj, file = proj_file);
}