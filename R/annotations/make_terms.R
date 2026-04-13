imports "annotation.workflow" from "seqtoolkit";

#' Generate Annotation Terms for Virtual Cell Modeling
#'
#' Executes the workflow step to parse DIAMOND BLASTP results and inject 
#' the extracted annotation terms into the virtual cell project session files. 
#' It supports two execution modes based on the "batch_process" configuration:
#' \itemize{
#'   \item \strong{Batch Mode}: Iterates through all available models, mapping 
#'     each model to its specific DIAMOND output directory to parse terms.
#'   \item \strong{Single Mode}: Processes a single project file, looking for 
#'     DIAMOND results in a shared workspace directory.
#' }
#'
#' @param app The application object or environment.
#' @param context The execution context for the current application run.
#'
#' @return Returns \code{invisible(NULL)}. Called for its side effects of 
#'   updating project session files with parsed annotation terms.
#'
#' @app make_terms
#' @export
[@app "make_terms"]
const make_terms = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));
    let diamond_workdir = WorkflowRender::workspace("make_diamond_hits");

    if (batch_process) {
        for(let model_dir in tqdm(list_batch_models())) {
            let proj_file = file.path(model_dir, "builder.gcproj");
            let blastp_dir = file.path(diamond_workdir, basename(model_dir));

            proj_file |> make_blastp_term(
                model_dir = blastp_dir
            );
        }
    } else {
        get_config("proj_file") |> make_blastp_term(
            model_dir = WorkflowRender::workspace("make_diamond_hits")
        );
    }

    invisible(NULL);
}

