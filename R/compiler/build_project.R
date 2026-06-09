#' Build GenBank Annotation Project into Virtual Cell Model
#'
#' Acts as the main entry point for building virtual cell models from
#' GenBank annotation project files. It automatically detects the processing
#' mode based on the \code{batch_process} configuration parameter:
#' \itemize{
#'   \item \strong{Batch Mode} (\code{batch_process = TRUE}): Iterates through
#'     all models returned by \code{\link{list_batch_models}}, reading the
#'     \code{builder.gcproj} file in each directory and saving the compiled
#'     output as \code{model.xml} (or a taxonomy-name-based XML file).
#'   \item \strong{Single Mode} (\code{batch_process = FALSE}): Processes a
#'     single project file specified in the configuration, allowing a specific
#'     virtual cell name (\code{vcell_name}) to be assigned.
#' }
#'
#' @param app The application object or environment, used for resolving
#'   workflow configuration via \code{get_config()}.
#' @param context The execution context object provided by the workflow engine.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item In batch mode: compiles each model into a GCMarkup XML file
#'       within its respective release subdirectory.
#'     \item In single mode: compiles the model into the XML file specified
#'       by the \code{model_file} configuration parameter.
#'   }
#'
#' @details
#' The function first loads the local data registry (required for constructing
#' the virtual cell component network), then dispatches to
#' \code{\link{compile_model}} for the actual compilation step.
#'
#' In batch mode, the output XML filename is derived from the organism's
#' taxonomy name (sanitized via \code{normalizeFileName}) when available,
#' falling back to the model directory basename otherwise.
#'
#' @seealso \code{\link{compile_model}} for the core compilation logic,
#'   \code{\link{list_batch_models}} for batch model directory listing.
#'
#' @examples
#' \dontrun{
#' # This function is typically invoked by the workflow engine:
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app build_project
#' @export
[@app "build_project"]
const build_project = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));
    let enzyme_fuzzy = as.logical(get_config("enzyme_fuzzy"));
    let registry = workflow::open_datapool(
        dir = get_config("localdb"), 
        enzyme_fuzzy = enzyme_fuzzy
    );

    registry <- set_kegg_pathways(registry, 
        maps = GCModeller::load_kegg_maps(), 
        reactions = GCModeller::kegg_reactions()
    );

    if (batch_process) {
        let gems_library_mode = as.logical(get_config("gems_library_mode"));
        let release_dir = get_config("gem_libout");

        for(let model_dir in list_batch_models()) {
            let proj_file = file.path(model_dir, "builder.gcproj");
            let save_xml = file.path(model_dir, "model.xml");

            if (gems_library_mode) {
                let tax_name = project::load(proj_file) |> scientific_name();
                let filename = normalizeFileName(tax_name, FALSE,maxchars = 64);

                save_xml <- file.path(release_dir, `${filename}.xml`);
            }

            message(`run build virtualcell model of '${basename(model_dir)}'!`);

            proj_file |> compile_model(
                save_model = save_xml, registry = registry);
        }

        message("finished run batch workflow for make virtual cell models!");
    } else {
        get_config("proj_file") |> compile_model(
            save_model = get_config("model_file"), 
            registry = registry, 
            # set vcell name only workspace in single genbank project mode
            vcell_name = get_config("vcell_name")   
        );
    }

    invisible(NULL);
}

