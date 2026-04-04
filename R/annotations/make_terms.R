imports "annotation.workflow" from "seqtoolkit";

#' make annotation terms of the cellular components inside the virtual cell model
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

#' write the diamond blastp result into the given virtualcell modelling project session file
#' 
#' @param model_dir the temp workspace of the virtual cell modelling
#' 
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