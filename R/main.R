#' Run the cell graph structure modelling workflow
#' 
#' @param src the ncbi genbank assembly source file its file path
#' @param outputdir the result outputdir and temp workspace
#'    location, default value is NULL means use the parent dir of
#'    the input src
#' @param biocyc the directory path to the biocyc reference pathway 
#' @param regprecise the data repository xml file path to the regprecise
#'    motif database.
#' 
#' @return this function has no return value, and then generated virtual
#'    cell model file will be saved at the ``outputdir`` with a fixed
#'    file name: ``model.vcell``. and the output dir also contains a html/pdf
#'    report about the virtual cell modelling result.
#'
const modelling_cellgraph = function(src, outputdir = NULL, 
                                     up_len = 150,
                                     biocyc = "./biocyc", 
                                     regprecise = "./RegPrecise.Xml") {
    WorkflowRender::init_context(outputdir || dirname(src));
    WorkflowRender::set_config(list(
        src        = normalizePath(src),
        biocyc     = normalizePath(biocyc),
        up_len     = up_len,
        regprecise = normalizePath(regprecise)
    ));
    WorkflowRender::run(registry = CellRender::annotation_workflow);
    WorkflowRender::finalize();

    invisible(NULL);
}

#' Run the kinetics parameter fitting based on the omics expression data
#' 
#' @param src the GCModeller virtual cell model assembly file its file path
#' 
const modelling_kinetics = function(src, outputdir = NULL) {
    invisible(NULL);
}