#' @title Run DIAMOND BlastP Search Workflow
#'
#' @description
#' Executes a DIAMOND blastp search workflow against multiple reference databases.
#' This function orchestrates the entire process: it configures the environment,
#' creates the necessary DIAMOND databases (if not already created), performs
#' the sequence alignment against EC numbers, subcellular locations, and
#' transcription factors, and manages the working directory state.
#'
#' The function supports two execution modes based on the \code{batch_process}
#' configuration:
#' \itemize{
#'   \item \strong{Batch Mode}: Iterates through all model subdirectories,
#'     locating each model's protein FASTA file and running DIAMOND searches
#'     in a model-specific workspace.
#'   \item \strong{Single Mode}: Processes a single protein FASTA file from
#'     the \code{make_genbank_proj} workspace output.
#' }
#'
#' @param app The application object passed by the workflow engine.
#'   Used for workflow context and configuration resolution via
#'   \code{get_config()} and \code{workfile()}.
#' @param context The execution context object provided by the workflow engine.
#'
#' @return This function is called for its side effects (creating .m8 output files).
#'   It does not return a value explicitly to the R session, but generates
#'   alignment result files in the DIAMOND workspace directory:
#'   \itemize{
#'     \item \code{ec_number.m8} - EC number annotation hits
#'     \item \code{subcellular.m8} - Subcellular localization hits
#'     \item \code{transcript_factor.m8} - Transcription factor hits
#'     \item \code{Pfam.csv} - Pfam domain annotation results
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Resolves the DIAMOND executable path from configuration or system PATH.
#'   \item Creates DIAMOND reference databases from the local data directory
#'     (EC numbers, subcellular locations, transcription factors) using
#'     \code{\link{make_diamond}}.
#'   \item Runs Pfam domain analysis via \code{\link{pfam_diamond}}.
#'   \item Executes DIAMOND BLASTP searches against each reference database.
#'   \item Restores the original working directory upon completion.
#' }
#'
#' @seealso \code{\link{make_diamond}} for database construction,
#'   \code{\link{pfam_diamond}} for Pfam domain analysis,
#'   \code{\link{make_terms}} for the downstream term-parsing step.
#'
#' @examples
#' \dontrun{
#' # This function is typically invoked by the workflow engine:
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app make_diamond_hits
#' @export
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
            "--threads", n_threads,
            "--outfmt","6", 
            "qseqid","stitle","pident","length","qstart","qend","sstart","send","evalue","bitscore"
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
        let source_dir = WorkflowRender::workspace("make_genbank_proj");
        let enable_blastp_cache = as.logical(get_config("enable_blastp_cache"));

        for(let model_dir in list_batch_models()) {
            let model_id = basename(model_dir);
            let proteins = file.path(source_dir, model_id, "proteins.fasta");
            let protein_pfam = file.path(dirname(proteins), "Pfam.csv");

            model_dir <- file.path(temp_dir, model_id);

            # create workspace dir for save diamond blastp result
            dir.create(model_dir);

            message(`make search for: ${proteins}`);
            message(`diamond blastp export to: ${model_dir}`);

            let ec_out = file.path(model_dir, "ec_number.m8");
            let cc_out = file.path(model_dir, "subcellular.m8");
            let tf_out = file.path(model_dir, "transcript_factor.m8");
            let check_size = 4[KB];
            
            if (file.size(protein_pfam) < check_size) {
                # no cache data
                # run process
                pfam_diamond(
                    proteins, 
                    workdir = dirname(proteins), 
                    diamond = diamond
                );
            }            

            if (enable_blastp_cache) {                
                let check_cache =  (file.size(ec_out) > check_size) 
                                && (file.size(cc_out) > check_size) 
                                && (file.size(tf_out) > check_size)
                ;

                if (check_cache) {
                    next;
                }
            }

            # then run diamond blastp search against the reference database
            diamond_blastp("ec_number", proteins, output = ec_out);
            diamond_blastp("subcellular", proteins, output = cc_out);
            diamond_blastp("transcript_factor",proteins, output = tf_out);

            message(`[${model_id}] diamond blastp search job done!`);
        }
    } else {
        # get genomics protein fasta sequence data file its file path for
        # run diamond blastp search
        let proteins = workfile("make_genbank_proj://proteins.fasta");
        let protein_pfam = file.path(dirname(proteins), "Pfam.csv");

        pfam_diamond(
            proteins, 
            workdir = dirname(proteins), 
            diamond = diamond
        );

        # then run diamond blastp search against the reference database
        diamond_blastp("ec_number", proteins, "ec_number.m8");
        diamond_blastp("subcellular", proteins, "subcellular.m8");
        diamond_blastp("transcript_factor",proteins, "transcript_factor.m8");
    }

    # restore the workdir finally.
    setwd(workdir);
}

