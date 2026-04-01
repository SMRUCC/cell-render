imports "bioseq.patterns" from "seqtoolkit";

[@app "make_TRN"]
const make_TRN = function(app, context) {
    if (check_build_module("TRN_network")) {
        let tfbs = workfile("tfbs_motif_scanning://tfbs_motifs.csv");
        let proj = project::load(get_config("proj_file"));

        tfbs = read.scans(tfbs);
        proj = proj |> set_tfbs(tfbs);

        project::save(proj, file = get_config("proj_file"));
    }
}