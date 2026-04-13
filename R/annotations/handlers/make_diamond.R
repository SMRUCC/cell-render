#' @title Build DIAMOND Reference Databases
#'
#' @description
#' Generates DIAMOND binary databases (`.dmnd`) from source FASTA files.
#' This utility function constructs three specific databases required for
#' the annotation workflow: Enzyme Commission (EC) numbers, Subcellular
#' locations, and Transcription Factors.
#'
#' @param local_db A character string specifying the directory path where
#'   the source FASTA files are located.
#' @param diamond A character string specifying the path to the DIAMOND
#'   executable. Defaults to \code{Sys.which("diamond")} if not provided.
#'
#' @return Invisible \code{NULL}. This function is called for its side effect
#'   of creating database files (e.g., \code{ec_number.dmnd}) in the current
#'   working directory.
#'
#' @details
#' This function expects the following source files to exist in \code{local_db}:
#' \itemize{
#'   \item \code{ec_numbers.fasta}
#'   \item \code{subcellular.fasta}
#'   \item \code{TF.fasta}
#' }
#' It utilizes \code{system2} to invoke the \code{diamond makedb} command.
#'
#' @importFrom base unlist
const make_diamond = function(local_db, diamond = Sys.which("diamond")) {
    let enzyme_db = file.path(local_db,"ec_numbers.fasta");
    let cc_location = file.path(local_db, "subcellular.fasta");
    let tf_db = file.path(local_db, "TF.fasta");

    diamond <- unlist(diamond);

    system2(diamond, c("makedb","--in",enzyme_db, "--db", "ec_number"), shell=TRUE);
    system2(diamond, c("makedb","--in",cc_location, "--db", "subcellular"), shell=TRUE);
    system2(diamond, c("makedb","--in",tf_db, "--db", "transcript_factor"), shell=TRUE);
}