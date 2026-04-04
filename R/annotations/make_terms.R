imports "annotation.workflow" from "seqtoolkit";

#' make annotation terms of the cellular components inside the virtual cell model
[@app "make_terms"]
const make_terms = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));

    if (batch_process) {
        for(let model_dir in tqdm(list_batch_models())) {
            let proj_file = file.path(model_dir, "builder.gcproj");
            let blastp_dir = file.path(WorkflowRender::workspace("make_diamond_hits"), basename(model_dir));

            proj_file |> make_blastp_term(
                model_dir = blastp_dir
            );
        }
    } else {
        get_config("proj_file") |> make_blastp_term(
            model_dir = WorkflowRender::workspace("make_diamond_hits")
        );
    }
}

const make_blastp_term = function(proj_file, model_dir) {
    # read the diamond blastp result files from the model worker dir
    let ec_number = read_m8(file.path(model_dir, "ec_number.m8"));
    let subcellular = read_m8(file.path(model_dir, "subcellular.m8"));
    let tf_list = read_m8(file.path(model_dir, "transcript_factor.m8"));
    let proj = project::load(proj_file);

    ec_number |> diamond_hitgroups |> set_blastp_result(proj, "ec_number");
    subcellular |> diamond_hitgroups |> set_blastp_result(proj, "subcellular_location");
    tf_list |> diamond_hitgroups |> set_blastp_result(proj, "transcript_factor");

    project::save(proj, file = proj_file);
}