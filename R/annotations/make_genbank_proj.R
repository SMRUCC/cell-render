imports "project" from "CellRender";
imports "workflow" from "CellRender";

#' Initialize a GenBank Project and Export Downstream Sequence Files
#'
#' @description
#' Loads GenBank source files specified in the application configuration,
#' initializes a new \code{project} object via \code{CellRender}, and saves the
#' resulting project state to disk. Additionally, it extracts and exports
#' essential FASTA files required for downstream bioinformatics analyses,
#' such as motif scanning and homology searches.
#'
#' This function serves as the workflow entry point and automatically detects
#' the processing mode from the \code{batch_process} configuration parameter:
#' \itemize{
#'   \item \strong{Batch Mode} (\code{batch_process = TRUE}): Scans the source
#'     directory for all GenBank files (\code{.gb}, \code{.gbk}, \code{.gbff})
#'     and processes each one individually via \code{\link{make_genbank_proj_file}}.
#'   \item \strong{Single Mode} (\code{batch_process = FALSE}): Processes a
#'     single GenBank source file specified by the \code{src} configuration key.
#' }
#'
#' @param app A CellRender workflow app object. Primarily used in this
#'   function to resolve working file paths via \code{workfile()}.
#' @param context A CellRender workflow context object.
#'   \emph{Note:} Currently unused in this specific function but required
#'   by the standard workflow execution interface.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item Creates a \code{builder.gcproj} project file in the release directory.
#'     \item Exports \code{upstream_locis.fasta} (TSS upstream sequences for motif
#'       scanning) and \code{proteins.fasta} (protein sequences for DIAMOND BLASTP)
#'       to the working directory.
#'   }
#'
#' @details
#' In batch mode, each GenBank file is processed in its own subdirectory
#' (named by the model accession ID) within both the release and working
#' directories. In single mode, files are written directly to the configured
#' paths without subdirectory nesting.
#'
#' @seealso \code{\link{make_genbank_proj_file}} for the core per-file processing
#'   logic, \code{\link{model_accession_id}} for accession ID extraction.
#'
#' @examples
#' \dontrun{
#' # This function is typically invoked by the workflow engine:
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app make_genbank_proj
#' @export
[@app "make_genbank_proj"]
const make_genbank_proj = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));
    let release_dir = get_config("release");
    let workdir = WorkflowRender::workspace(app);

    if (batch_process) {  
        let genbank_files = list.files(get_config("src"), 
                                pattern = c("*.gb","*.gbk","*.gbff"), recursive = TRUE);

        message(`Build virtualcell community model based ${length(genbank_files)} genbank source files!`);

        for(let file in tqdm(genbank_files)) {
            file |> make_genbank_proj_file(
                release_dir = release_dir,
                workdir = workdir,
                batch_process = TRUE
            );
        }
    } else {
        # input source is a single genbank file
        let gb_src = get_config("src");        

        gb_src |> make_genbank_proj_file(
            release_dir = release_dir,
            workdir = workdir,
            batch_process = FALSE
        );
    }

    invisible(NULL);
}

