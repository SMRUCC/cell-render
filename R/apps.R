const workflow_registry = function() {
    WorkflowRender::hook(app("extract_genbank_src", extract_gene_table));
} 