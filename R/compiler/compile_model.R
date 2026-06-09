#' Compile Project File into Virtual Cell Model XML
#'
#' Core processing function that loads a GenBank annotation project, constructs
#' the virtual cell component network using a local data registry, and exports
#' the resulting model as a GCMarkup XML file.
#'
#' @param proj_file Character. The file path to the GenBank annotation
#'   project file (e.g., \code{.gcproj}).
#' @param save_model Character. The destination file path where the generated
#'   virtual cell model will be saved. The model is saved in GCMarkup XML format.
#' @param registry A local data repository (data pool object) containing the
#'   necessary reference data required to construct the virtual cell component
#'   network.
#' @param vcell_name Character, optional. The specific name to assign to the
#'   generated virtual cell. If \code{NULL} (default), the naming defaults to
#'   the project's internal settings. This parameter is typically utilized only
#'   in single-project mode.
#'
#' @return Writes an XML file to the specified \code{save_model} path and
#'   returns \code{invisible(NULL)}.
#'
#' @details
#' The compilation process:
#' \enumerate{
#'   \item Loads the project file from \code{proj_file}.
#'   \item Constructs the virtual cell model via \code{build()}, using the
#'     provided \code{registry} as a reference data pool and optionally
#'     setting the virtual cell name.
#'   \item Serializes the virtual cell model to GCMarkup XML format.
#'   \item Writes the XML content to the file specified by \code{save_model}.
#' }
#'
#' @seealso \code{\link{build_project}} which calls this function,
#'   \code{\link{make_terms}} and \code{\link{make_TRN}} which produce the
#'   annotation data consumed during compilation.
#'
#' @examples
#' \dontrun{
#' # Compile a single project into a virtual cell model
#' compile_model(
#'   proj_file = "release/builder.gcproj",
#'   save_model = "release/model.xml",
#'   registry = CellRender::load_registry(),
#'   vcell_name = "E_coli_K12"
#' )
#' }
#'
#' @keywords internal
#' @export
const compile_model = function(proj_file, save_model, registry, vcell_name = NULL) {
    let proj = project::load(proj_file);
    let vcell = proj |> build(
        datapool   = registry, 
        vcell_name = vcell_name
    );

    # write the virtual cell model into a xml file
    vcell 
    |> xml 
    |> writeLines(con = save_model)
    ; 
}