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
#' @param proj_file Character. The file path to the GenBank annotation 
#'   project file (\code{.gcproj}) to be updated.
#' @param model_dir Character. The file path to the temporary workspace 
#'   directory containing the DIAMOND BLASTP \code{.m8} result files.
#'
#' @return Returns \code{invisible(NULL)}. The project file specified by 
#'   \code{proj_file} is modified and saved with the new annotation terms.
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

    project::save(proj, file = proj_file);
}