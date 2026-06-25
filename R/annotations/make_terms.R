imports "annotation.workflow" from "seqtoolkit";

#' Generate Annotation Terms for Virtual Cell Modeling
#'
#' Executes the workflow step to parse DIAMOND BLASTP results and inject
#' the extracted annotation terms into the virtual cell project session files.
#' It supports two execution modes based on the \code{batch_process} configuration:
#' \itemize{
#'   \item \strong{Batch Mode}: Iterates through all available models, mapping
#'     each model to its specific DIAMOND output directory to parse terms.
#'   \item \strong{Single Mode}: Processes a single project file, looking for
#'     DIAMOND results in a shared workspace directory.
#' }
#'
#' @param app The application object or environment, used for resolving
#'   workflow configuration via \code{get_config()}.
#' @param context The execution context for the current application run,
#'   provided by the \pkg{WorkflowRender} engine.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item Parses \code{ec_number.m8}, \code{subcellular.m8}, and
#'       \code{transcript_factor.m8} files from the DIAMOND output directory.
#'     \item Writes grouped annotation terms into the project session file
#'       (\code{builder.gcproj}).
#'     \item Exports an enzyme annotation table as \code{enzymes.csv} in the
#'       project directory.
#'   }
#'
#' @details
#' The parsed annotation terms serve three distinct cellular network
#' construction purposes:
#' \itemize{
#'   \item \code{ec_number} terms feed into the \strong{metabolic network}.
#'   \item \code{subcellular_location} terms feed into the \strong{transmembrane network}.
#'   \item \code{transcript_factor} terms feed into the \strong{gene expression
#'     transcription regulation network}.
#' }
#'
#' @seealso \code{\link{make_blastp_term}} for the core term-parsing logic,
#'   \code{\link{make_diamond_hits}} for the preceding DIAMOND search step,
#'   \code{\link{list_batch_models}} for batch model directory listing.
#'
#' @examples
#' \dontrun{
#' # This function is typically invoked by the workflow engine:
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app make_terms
#' @export
[@app "make_terms"]
const make_terms = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));
    let diamond_workdir = WorkflowRender::workspace("make_diamond_hits");
    let release_dir = get_config("release");

    if (batch_process) {
        for(let model_dir in tqdm(list_batch_models())) {
            let proj_file = file.path(model_dir, "builder.gcproj");
            let model_id = basename(model_dir);
            let blastp_dir = file.path(diamond_workdir, model_id );
            let traits = metaTraits(workfile(`make_genbank_proj://${model_id}/proteins.tsv`));

            proj_file |> make_blastp_term(
                traits = traits,
                model_dir = blastp_dir
            );

            write.csv(as.data.frame(traits ), file = file.path(release_dir, model_id, "metaTraits.csv"));
        }
    } else {
        get_config("proj_file") |> make_blastp_term(
            model_dir = WorkflowRender::workspace("make_diamond_hits")
        );
    }

    invisible(NULL);
}

