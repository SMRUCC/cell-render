#' Run the modelling workflow
#' 
#' @param src the ncbi genbank assembly source file its file path
#' @param outputdir the result outputdir and temp workspace
#'    location, default value is NULL means use the parent dir of
#'    the input src
#' 
const run_workflow = function(src, outputdir = NULL) {
    WorkflowRender::init_context(outputdir || dirname(src));
    WorkflowRender::run(registry = CellRender::workflow_registry);
    WorkflowRender::finalize();
}