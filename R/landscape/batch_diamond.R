#' Batch DIAMOND BLASTP Search Against Enzyme Database
#'
#' A helper function that formats a local EC number reference database and
#' executes DIAMOND BLASTP in batch mode for all protein FASTA files in a
#' given source directory.
#'
#' @param source_dir A character string specifying the directory path containing
#'   input protein FASTA files (\code{.faa}).
#' @param result_dir A character string specifying the directory path where the
#'   DIAMOND output files (\code{.m8} format) will be saved.
#' @param diamond A character string specifying the file path to the DIAMOND
#'   executable. Defaults to \code{Sys.which("diamond")}.
#' @param n_threads An integer specifying the number of CPU threads to allocate
#'   to DIAMOND. Defaults to 32.
#'
#' @details
#' This function temporarily changes the working directory to \code{result_dir}
#' for the DIAMOND database creation and search operations, then restores the
#' original working directory upon completion.
#'
#' The EC number reference database (\code{ec_numbers.fasta}) is sourced from
#' the package's built-in data directory. A DIAMOND database is created
#' on-the-fly before running the BLASTP searches.
#'
#' For each \code{.faa} file found in \code{source_dir}, a corresponding
#' \code{.m8} alignment result file is generated in \code{result_dir}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of creating DIAMOND database files and BLASTP alignment output files
#'   (\code{.m8}) in the \code{result_dir}.
#'
#' @seealso \code{\link{make_diamond_hits}} for the main DIAMOND workflow step,
#'   \code{\link{make_diamond}} for building individual reference databases.
#'
#' @examples
#' \dontrun{
#' # Run batch DIAMOND BLASTP on all .faa files
#' batch_diamond(
#'   source_dir = "data/proteins",
#'   result_dir = "results/diamond",
#'   diamond = "/usr/bin/diamond",
#'   n_threads = 16
#' )
#' }
#'
#' @keywords internal
#' @export
const batch_diamond = function(source_dir, result_dir, 
                               diamond = Sys.which("diamond"), 
                               n_threads = 32) {
    let current_dir = getwd();
    # use the ec_numbers.fasta as reference db
    # which is inside the data dir of current R package
    let local_db = file.path(@datadir, "ec_numbers.fasta");

    setwd(blastp_dir);        
    system2(diamond, c("makedb","--in", local_db, "--db", "ec_number"), shell=TRUE);

    for(let faa in list.files(source_dir, pattern = "*.faa")) {
        system2(diamond, c("blastp",
            "--db", "ec_number", 
            "--query", faa, 
            "--out", file.path(result_dir, `${basename(faa)}.m8`), 
            "--threads", n_threads
        ), shell=TRUE)
        ;
    }

    setwd(current_dir);
}