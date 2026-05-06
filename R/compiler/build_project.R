#' Build GenBank Annotation Project into Virtual Cell Model
#'
#' Acts as the main entry point for building virtual cell models from
#' GenBank annotation project files. It automatically detects the processing
#' mode based on the "batch_process" configuration parameter:
#' \itemize{
#'   \item \strong{Batch Mode} (\code{batch_process = TRUE}): Iterates through 
#'     all models returned by \code{list_batch_models()}, reading the 
#'     \code{builder.gcproj} file in each directory and saving the compiled 
#'     output as \code{model.xml}.
#'   \item \strong{Single Mode} (\code{batch_process = FALSE}): Processes a 
#'     single project file specified in the configuration, allowing a specific
#'     virtual cell name (\code{vcell_name}) to be assigned.
#' }
#'
#' @param app The application object or environment.
#' @param context The execution context for the current application run.
#'
#' @return Returns \code{invisible(NULL)}. This function is called primarily 
#'   for its side effects: generating virtual cell model XML files and logging 
#'   progress messages to the console.
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
        maps = GCModeller::kegg_maps(), 
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

