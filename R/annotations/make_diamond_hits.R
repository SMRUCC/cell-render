#' @title Run DIAMOND BlastP Search Workflow
#'
#' @description
#' Executes a DIAMOND blastp search workflow against multiple reference databases.
#' This function orchestrates the entire process: it configures the environment,
#' creates the necessary DIAMOND databases (if not already created), performs
#' the sequence alignment against EC numbers, subcellular locations, and
#' transcription factors, and manages the working directory state.
#'
#' @param app The application object passed by the workflow engine. 
#'   Used for workflow context and configuration.
#' @param context The execution context object provided by the workflow engine.
#'
#' @return This function is called for its side effects (creating .m8 output files).
#'   It does not return a value explicitly to the R session, but generates
#'   alignment result files in the workflow workspace.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Retrieves configuration for the DIAMOND path, database directory,
#'         and thread count.
#'   \item Sets the working directory to the specific workflow workspace.
#'   \item Calls \code{make_diamond} to build binary databases from FASTA sources.
#'   \item Runs DIAMOND blastp against three databases: 
#'         \itemize{
#'           \item \strong{ec_number}: Enzyme Commission numbers.
#'           \item \strong{subcellular}: Subcellular location predictions.
#'           \item \strong{transcript_factor}: Transcription factors.
#'         }
#'   \item Restores the original working directory upon completion.
#' }
#'
#' @seealso \code{\link{make_diamond}} for the database creation logic.
#'
[@app "make_diamond_hits"]
const make_diamond_hits = function(app, context) {
    let diamond = get_config("diamond");  # diamond path
    let localdb = get_config("localdb");  # reference database dir path
    
    # run this diamond blastp alignment workflow in batch process mode?
    let batch_process = as.logical(get_config("batch_process"));
    let n_threads = get_config("n_threads");
    let workdir = getwd();

    # a helper wrapper function of the diamond blastp search 
    # commandline calls
    let diamond_blastp = function(db, proteins, output) {
        system2(diamond, c("blastp",
            "--db", db, 
            "--query", proteins, 
            "--out", output, 
            "--threads", n_threads
        ), shell=TRUE)
        ;
    }
    let temp_dir = WorkflowRender::workspace(app);

    # set current workdir to the temp workspace of 
    # this `make_diamond_hits` workflow module
    setwd(temp_dir);
    # make reference database
    make_diamond(localdb, diamond);

    if (batch_process) {
        for(let model_dir in list.dirs(WorkflowRender::workspace("make_genbank_proj"), recursive = FALSE)) {
            let proteins = file.path(model_dir, "proteins.fasta");
            let model_id = basename(model_dir);

            model_dir <- file.path(temp_dir, model_id);

            # then run diamond blastp search against the reference database
            diamond_blastp("ec_number", proteins, file.path(model_dir, "ec_number.m8"));
            diamond_blastp("subcellular", proteins, file.path(model_dir, "subcellular.m8"));
            diamond_blastp("transcript_factor",proteins, file.path(model_dir, "transcript_factor.m8"));
        }
    } else {
        # get genomics protein fasta sequence data file its file path for
        # run diamond blastp search
        let proteins = workfile("make_genbank_proj://proteins.fasta");

        # then run diamond blastp search against the reference database
        diamond_blastp("ec_number", proteins, "ec_number.m8");
        diamond_blastp("subcellular", proteins, "subcellular.m8");
        diamond_blastp("transcript_factor",proteins, "transcript_factor.m8");
    }

    # restore the workdir finally.
    setwd(workdir);
}

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