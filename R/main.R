#' Run the cell graph structure modelling workflow
#' 
#' @param src the ncbi genbank assembly source file its file path
#' @param outputdir the result outputdir and temp workspace
#'    location, default value is NULL means use the parent dir of
#'    the input src
#' @param biocyc the directory path to the local reference dataset for make the cellular component annotation. 
#' @param up_len an integer vector for specific the length of the TSS upstream site
#' 
#' @return this function has no return value, and then generated virtual
#'    cell model file will be saved at the ``outputdir`` with a fixed
#'    file name: ``model.vcell``. and the output dir also contains a html/pdf
#'    report about the virtual cell modelling result.
#'
const modelling_cellgraph = function(src, outputdir = NULL, 
                                     up_len = 150, 
                                     localdb = NULL, 
                                     diamond = Sys.which("diamond"), 
                                     domain = c("bacteria", "plant", "animal", "fungi"),
                                     builds = c("TRN_network","Metabolic_network"),
                                     n_threads = 32) {

    WorkflowRender::init_context(outputdir || dirname(src));
    WorkflowRender::set_config(list(
        src        = normalizePath(src),
        localdb    = localdb || normalizePath(@datadir),
        up_len     = up_len,
        diamond    = unlist(diamond),
        n_threads  = n_threads,
        domain     = .Internal::first(domain),
        builds     = builds,
        release    = file.path(workdir_root(), "release")
        proj_file  = file.path(workdir_root(), "release", "builder.gcproj"),
        model_file = file.path(workdir_root(), "release", "model.xml") 
    ));
    WorkflowRender::run(registry = CellRender::annotation_workflow);
    WorkflowRender::finalize();

    invisible(NULL);
}
