#' List All Batch Model Directories
#'
#' Retrieves the list of model subdirectories from the release directory
#' configured in the current workflow session. Each subdirectory represents
#' a separate virtual cell model that was created during batch processing.
#'
#' @return A character vector of directory paths, one per model. Returns an
#'   empty vector if no model directories exist.
#'
#' @details
#' This function reads the \code{release} configuration parameter to
#' determine the root directory, then lists all immediate subdirectories.
#' Each subdirectory is expected to contain a \code{builder.gcproj} file
#' produced by \code{\link{make_genbank_proj_file}}.
#'
#' @seealso \code{\link{make_genbank_proj}} for the workflow step that
#'   creates the model directories, \code{\link{build_project}} and
#'   \code{\link{make_terms}} which iterate over the returned directories.
#'
#' @examples
#' \dontrun{
#' # List all models in the release directory
#' models <- list_batch_models()
#' print(models)
#' # [1] "/path/to/release/GCF_000005845" "/path/to/release/GCF_000006945"
#' }
#'
#' @keywords internal
#' @export
const list_batch_models = function() {
    list.dirs(get_config("release"), recursive = FALSE);
}