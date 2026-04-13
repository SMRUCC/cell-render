#' Register the virtual cell network modelling workflows
#' 
#' @return this function has no return value
#' 
const annotation_workflow = function() {
    print("Config the cellular graph network annotation & modelling workflow...");

    # steps for make genbank annotation project
    WorkflowRender::hook(CellRender::make_genbank_proj);
    WorkflowRender::hook(CellRender::tfbs_motif_scanning);
    # blastp search and then make annotation terms
    WorkflowRender::hook(CellRender::make_diamond_hits);
    WorkflowRender::hook(CellRender::make_terms);
    WorkflowRender::hook(CellRender::make_TRN);

    # steps for create virtual cell model from genbank annotation project file
    WorkflowRender::hook(CellRender::build_project);

    WorkflowRender::summary();

    invisible(NULL);
} 
