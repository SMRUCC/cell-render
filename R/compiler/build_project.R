#' build the genbank annotation project file into the virtual cell model
[@app "build_project"]
const build_project = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));
    let registry = open_datapool(dir = get_config("localdb"));

    if (batch_process) {
        for(let model_dir in list_batch_models()) {
            let proj_file = file.path(model_dir, "builder.gcproj");
            let save_xml = file.path(model_dir, "model.xml");

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

#' compile the project file into virtual cell model file
#' 
#' @param proj_file the file path to the genbank annotation project file
#' @param save_model the file save path of the generated virtual cell model file, model file will be saved in GCMarkup xml file format.
#' @param registry a local data repository that contains the necessary data for make the virtual cell component network.
#' 
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