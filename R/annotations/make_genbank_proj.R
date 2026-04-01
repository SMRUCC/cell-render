imports "project" from "CellRender";
imports "workflow" from "CellRender";

#' Extract the genbank source
#' 
#' @param app the current workflow app object
#' @param context the workflow context object
#' 
[@app "make_genbank_proj"]
const make_genbank_proj = function(app, context) {
    let gb_src = load_genbanks(get_config("src"));
    let proj = project::new(gb_src);

    project::save(proj, file = get_config("proj_file"));

    write.fasta(tss_upstream(proj), file = workfile(app, "upstream_locis.fasta")); 
    workflow::save_proteins(proj, file = workfile(app, "proteins.fasta"));
}