imports "bioseq.patterns" from "seqtoolkit";

#' Build the Transcription Regulation Network (TRN)
#'
#' Constructs the transcription regulation network within the virtual cell
#' model by integrating the TFBS motif scanning results into the project
#' session. This step is conditionally executed and will only run if the
#' \code{TRN_network} module is enabled in the virtual cell build configuration.
#'
#' The TRN is built by combining two data sources:
#' \enumerate{
#'   \item \strong{TFBS motif sites} — obtained from the preceding
#'     \code{\link{tfbs_motif_scanning}} workflow step, which identifies
#'     potential transcription factor binding sites in gene upstream regions.
#'   \item \strong{Transcription factor annotations} — obtained from the
#'     \code{\link{make_terms}} workflow step, which identifies which proteins
#'     are transcription factors via DIAMOND BLASTP.
#' }
#'
#' @param app The application object passed by the workflow engine.
#'   Used for resolving workflow file paths via \code{workfile()}.
#' @param context The execution context object provided by the workflow engine.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item Reads the TFBS motif scan results from the workspace.
#'     \item Injects the TFBS data into the project session via \code{set_tfbs}.
#'     \item Saves the updated project file to disk.
#'   }
#'
#' @details
#' This function checks whether the \code{TRN_network} build module is enabled
#' via \code{\link{check_build_module}}. If not enabled, the function exits
#' immediately without performing any operations.
#'
#' @seealso \code{\link{tfbs_motif_scanning}} for the preceding TFBS scanning step,
#'   \code{\link{make_terms}} for the transcription factor annotation step,
#'   \code{\link{check_build_module}} for module enablement checking.
#'
#' @examples
#' \dontrun{
#' # This function is typically invoked by the workflow engine:
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app make_TRN
#' @export
[@app "make_TRN"]
const make_TRN = function(app, context) {
    if (check_build_module("TRN_network")) {
        let tfbs = workfile("tfbs_motif_scanning://tfbs_motifs.csv");
        let proj = project::load(get_config("proj_file"));

        tfbs = read.scans(tfbs);
        proj = proj |> set_tfbs(tfbs);

        project::save(proj, file = get_config("proj_file"));
    }
}