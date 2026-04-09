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

#' Compile Project File into Virtual Cell Model XML
#'
#' Core processing function that loads a GenBank annotation project, constructs 
#' the virtual cell component network using a local data registry, and exports 
#' the resulting model as a GCMarkup XML file.
#'
#' @param proj_file Character. The file path to the GenBank annotation 
#'   project file (e.g., \code{.gcproj}).
#' @param save_model Character. The destination file path where the generated 
#'   virtual cell model will be saved. The model is saved in GCMarkup XML format.
#' @param registry A local data repository (data pool object) containing the 
#'   necessary reference data required to construct the virtual cell component 
#'   network.
#' @param vcell_name Character, optional. The specific name to assign to the 
#'   generated virtual cell. If \code{NULL} (default), the naming defaults to 
#'   the project's internal settings. This parameter is typically utilized only 
#'   in single-project mode.
#'
#' @return Writes an XML file to the specified \code{save_model} path and 
#'   returns \code{invisible(NULL)}.
#'
#' @keywords internal
#' @export
const compile_model = function(proj_file, save_model, registry, vcell_name = NULL) {
    let proj = project::load(proj_file);
    let vcell = proj |> build(
        datapool   = registry, 
        vcell_name = vcell_name
    );

    # write the virtual cell model into a xml file
    vcell 
    |> xml 
    |> writeLines(con = save_model)
    ; 
}