#' Register the virtual cell network modelling workflows
#' 
#' @return this function has no return value
#' 
const annotation_workflow = function() {
    print("Config the cellular graph network annotation & modelling workflow...");

    # steps for make genbank annotation project
    WorkflowRender::hook(make_genbank_proj);
    WorkflowRender::hook(tfbs_motif_scanning);
    WorkflowRender::hook(make_diamond_hits);
    WorkflowRender::hook(make_terms);

    # steps for create virtual cell model from genbank annotation project file
    WorkflowRender::hook(build_project);

    invisible(NULL);
} 