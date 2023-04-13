#' Register the virtual cell network modelling workflows
#' 
const annotation_workflow = function() {
    WorkflowRender::hook(app("extract_genbank_src", extract_gene_table));

    invisible(NULL);
} 