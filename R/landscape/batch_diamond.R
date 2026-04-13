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
#' to generate the DIAMOND database (\code{ec_number.dmnd}) using an internal 
#' package reference file (\code{ec_numbers.fasta}). It then iterates through 
#' all \code{.faa} files in \code{source_dir} and executes a BLASTP search. 
#' The original working directory is restored upon completion.
#'
#' @return Returns \code{invisible(NULL)}. Output \code{.m8} files and the 
#'   DIAMOND database file are written to \code{result_dir}.
#'
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