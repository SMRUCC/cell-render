#' Register the virtual cell network modelling workflows
#' 
const annotation_workflow = function() {
    WorkflowRender::hook(app("extract_genbank_src", extract_gene_table));
    WorkflowRender::hook(app("extract_tfbs_motifs", tfbs_motif_scanning));
    WorkflowRender::hook(app("assemble_TRN", assemble_transcript_graph));
    WorkflowRender::hook(app("assemble_RXN", assemble_metabolic_graph));

    invisible(NULL);
} 