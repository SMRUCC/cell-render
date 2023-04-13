#' Run the cell graph structure modelling workflow
#' 
#' @param src the ncbi genbank assembly source file its file path
#' @param outputdir the result outputdir and temp workspace
#'    location, default value is NULL means use the parent dir of
#'    the input src
#' @param biocyc the directory path to the biocyc reference pathway 
#' 
const modelling_cellgraph = function(src, outputdir = NULL, biocyc = "./biocyc") {
    WorkflowRender::init_context(outputdir || dirname(src));
    WorkflowRender::set_config(list(
        src = normalizePath(src)
    ));
    WorkflowRender::run(registry = CellRender::annotation_workflow);
    WorkflowRender::finalize();
}

#' Run the kinetics parameter fitting based on the omics expression data
#' 
#' @param src the GCModeller virtual cell model assembly file its file path
#' 
const modelling_kinetics = function(src, outputdir = NULL) {

}