#' Register the virtual cell network modelling workflows
#' 
#' @return this function has no return value
#' 
const annotation_workflow = function() {
    print("Config the cellular graph network annotation & modelling workflow...");

    WorkflowRender::hook(make_genbank_proj);
    WorkflowRender::hook(tfbs_motif_scanning);
    WorkflowRender::hook(make_diamond_hits);
    WorkflowRender::hook(app("assemble_RXN", assemble_metabolic_graph));

    invisible(NULL);
} 