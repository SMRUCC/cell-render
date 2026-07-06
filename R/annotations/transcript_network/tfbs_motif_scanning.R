imports "bioseq.fasta" from "seqtoolkit";
imports "bioseq.patterns" from "seqtoolkit";

#' Scan for Transcription Factor Binding Sites (TFBS) Motifs
#' 
#' @description Performs motif scanning on gene upstream sequences to identify
#'   Transcription Factor Binding Sites (TFBS). This function is conditionally
#'   executed and will only run if the \code{TRN_network} module is enabled in
#'   the virtual cell build configuration.
#'
#' @param app The application environment or object, used to determine the
#'   output working directory via \code{workfile()}.
#' @param context The execution context passed by the GCModeller framework.
#'
#' @details The function executes a specific workflow for regulatory network
#'   construction:
#'   \enumerate{
#'     \item \strong{Input Retrieval}: Reads the TSS upstream sequence FASTA
#'       file produced by \code{\link{make_genbank_proj}}.
#'     \item \strong{Module Check}: Verifies that the \code{TRN_network} build
#'       module is enabled via \code{\link{check_build_module}}. If not, the
#'       function exits early.
#'     \item \strong{Motif Scanning}: Loads the TFBS motif database from the
#'       package data directory and runs \code{GCModeller::scan_motifs} with
#'       a 90\% identity cutoff.
#'     \item \strong{Result Export}: Writes the motif scan results as a CSV file
#'       (\code{tfbs_motifs.csv}) for downstream TRN construction by
#'       \code{\link{make_TRN}}.
#'   }
#'
#'   The output CSV contains motif hit records with the following structure:
#'   \itemize{
#'     \item \code{motif.id} — The motif identifier from the reference database.
#'     \item \code{family} — The transcription factor family classification.
#'     \item \code{name} — The descriptive name of the motif/TF.
#'     \item \code{TF.name} — The gene locus tag of the matched transcription factor.
#'     \item \code{protein_model} — The protein model identifier.
#'   }
#'
#' @return \code{invisible(NULL)}. This function is called for its side effect
#'   of writing the \code{tfbs_motifs.csv} file to the workflow workspace.
#'
#' @seealso \code{\link{make_genbank_proj}} for upstream sequence extraction,
#'   \code{\link{make_TRN}} for the downstream TRN construction step,
#'   \code{\link{check_build_module}} for module enablement checking.
#'
#' @examples
#' \dontrun{
#' # This function is typically invoked by the workflow engine:
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app tfbs_motif_scanning
#' @export
[@app "tfbs_motif_scanning"]
const tfbs_motif_scanning = function(app, context) {    
    # run this diamond blastp alignment workflow in batch process mode?
    let batch_process = as.logical(get_config("batch_process"));
    let n_threads = get_config("n_threads");
    # get the motif database directory that contains multiple meme motif model files
    let motifs_db = file.path(get_config("localdb"), "motifs", get_config("domain"));

    # run motif site scanning if the TRN_network module build in virtual cell is enabled
    if (check_build_module("TRN_network")) {
        if (batch_process) {
            let source_dir = WorkflowRender::workspace("make_genbank_proj");

            for(let model_dir in list_batch_models()) {
                # set the TSS upstream region site fasta file
                let upstream_seq = file.path(model_dir, "upstream_locis.fasta");
                let outfile = file.path(model_dir, "tfbs_motifs.csv");
                # make TFBS site scanning on the TSS upstream region sites
                # search site against the reference motif search.
                let motifs = GCModeller::scan_motifs(
                    db = motifs_db, 
                    seqs = upstream_seq, 
                    workdir = model_dir, 
                    n_threads = n_threads,
                    identities_cutoff = 0.9
                );

                write.csv(motifs, file = outfile, row.names = FALSE);
            }
        } else {
            # set the TSS upstream region site fasta file
            let upstream_seq = workfile("make_genbank_proj://upstream_locis.fasta");
            
            # make TFBS site scanning on the TSS upstream region sites
            # search site against the reference motif search.
            let motifs = GCModeller::scan_motifs(
                db = motifs_db, 
                seqs = upstream_seq, 
                workdir = workfile(app, "tfbs_sites"), 
                n_threads = get_config("n_threads"),
                identities_cutoff = 0.9
            );
            let outfile = workfile(app, "tfbs_motifs.csv");

            # seeds 
            # $"{motif.id} {motif.family} [{motif.name}]"
            #
            #         motif.id family [name]
            #              --- ---- -------------------   
            # "AsnC; MOTIF 557 AsnC [YwrC - Bacillales] Motif_1" 

            # TF 
            #
            #    TF.name motif.id protein_model
            #     ------ ----     --------------
            #   > jk0144 1269 Protein:144-jk0144

            write.csv(motifs, file = outfile, row.names = FALSE);
        }
    }
}

