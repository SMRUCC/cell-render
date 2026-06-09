require(GCModeller);
require(WorkflowRender);

#' Package Initialization Hook
#'
#' Executed automatically when the CellRender package is loaded via
#' \code{library(CellRender)}. It ensures that the core dependency packages
#' \pkg{GCModeller} and \pkg{WorkflowRender} are available in the current
#' R session, as they provide the foundational workflow engine and rendering
#' infrastructure required by all virtual cell modelling pipelines.
#'
#' @param libname Character. The library directory from which the package is
#'   being loaded. Supplied by the R runtime.
#' @param pkgname Character. The name of the package being loaded. Supplied
#'   by the R runtime.
#'
#' @return This function returns \code{invisible(NULL)}. It is called
#'   exclusively for its side effect of loading dependency packages.
#'
#' @keywords internal
const .onLoad = function() {
    
}