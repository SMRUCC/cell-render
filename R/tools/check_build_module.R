#' Check Whether a Build Module is Enabled
#'
#' Checks if a specific virtual cell build module is enabled in the
#' current workflow configuration. This is used to conditionally execute
#' workflow steps based on user-specified build targets.
#'
#' @param flag Character. The name of the build module to check
#'   (e.g., \code{"TRN_network"}, \code{"metabolic_network"}).
#'   The comparison is case-insensitive.
#'
#' @return Logical. \code{TRUE} if the specified module is present in the
#'   \code{builds} configuration, \code{FALSE} otherwise.
#'
#' @details
#' The function reads the \code{builds} configuration parameter (a character
#' vector of enabled module names) and performs a case-insensitive match
#' against the provided \code{flag}.
#'
#' @seealso \code{\link{make_TRN}} and \code{\link{tfbs_motif_scanning}}
#'   which use this function to conditionally execute the TRN network
#'   construction steps.
#'
#' @examples
#' \dontrun{
#' # Check if the TRN network module should be built
#' if (check_build_module("TRN_network")) {
#'   message("TRN network module is enabled")
#' }
#' }
#'
#' @keywords internal
#' @export
const check_build_module = function(flag) {
    let builds = get_config("builds");
    let check = any(tolower(flag) == tolower(builds));

    return(check);
}