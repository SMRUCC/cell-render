#' Extract the genbank source
#' 
#' @param app the current workflow app object
#' @param context the workflow context object
#' 
const extract_gene_table = function(app, context) {
    const src = context$configs$src;    
    # get workspace dir path of current workflow app module
    const workdir = WorkflowRender::workspace(app);
    const verbose = as.logical(getOption("verbose"));
    const upstream_size = context$configs$up_len || 150;

    if (verbose) {
        print("app workspace for extract gene table:");
        print(workdir);
    }

    CellRender::extract_gbff(src, workdir, upstream_size, 
        verbose = verbose);
}