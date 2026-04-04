const list_batch_models = function() {
    list.dirs(WorkflowRender::workspace("make_genbank_proj"), recursive = FALSE);
}