#' build project into the virtual cell model
[@app "build_project"]
const build_project = function(app, context) {
    let proj = project::load(get_config("proj_file"));
    let registry = open_datapool(dir = get_config("localdb"));
    let vcell_name = get_config("vcell_name");
    let vcell = proj |> build(
        datapool   = registry, 
        vcell_name = vcell_name
    );

    vcell 
    |> xml 
    |> writeLines(con = get_config("model_file"))
    ; 
}