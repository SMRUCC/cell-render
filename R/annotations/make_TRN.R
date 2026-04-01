imports "bioseq.patterns" from "seqtoolkit";

[@app "make_TRN"]
const make_TRN = function(app, context) {
    let tfbs = workfile("tfbs_motif_scanning://tfbs_motifs.csv");
    let proj = project::load(workfile("make_terms://builder.gcproj"));

    tfbs = read.scans(tfbs);
    proj = proj |> set_tfbs(tfbs);

    project::save(proj, file = workfile(app, "builder.gcproj"));
}