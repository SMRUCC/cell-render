imports "annotation.workflow" from "seqtoolkit";

[@app "make_terms"]
const make_terms = function(app, context) {
    let ec_number = read_m8(workfile("make_diamond_hits://ec_number.m8"));
    let subcellular = read_m8(workfile("make_diamond_hits://subcellular.m8"));
    let tf_list = read_m8(workfile("make_diamond_hits://transcript_factor.m8"));
    let proj = project::load(workfile("make_genbank_proj://builder.gcproj"));

    ec_number |> diamond_hitgroups |> set_blastp_result(proj, "ec_number");
    subcellular |> diamond_hitgroups |> set_blastp_result(proj, "subcellular_location");
    tf_list |> diamond_hitgroups |> set_blastp_result(proj, "transcript_factor");

    project::save(proj, file = workfile(app, "builder.gcproj"));
}