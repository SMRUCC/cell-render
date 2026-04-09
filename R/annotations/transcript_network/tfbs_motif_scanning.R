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
#'     \item \strong{Input Retrieval}: Reads the TSS upstream region sequences 
#'       (generated from \code{make_genbank_proj}) as a FASTA file.
#'     \item \strong{Database Query}: Constructs the path to the local MEME 
#'       motif database using the \code{localdb} and \code{domain} configurations.
#'     \item \strong{Scanning}: Calls \code{GCModeller::scan_motifs} to scan 
#'       the upstream sequences against the motif database. The search uses 
#'       a strict identity cutoff of \code{0.9} and utilizes multithreading 
#'       based on the \code{n_threads} configuration.
#'     \item \strong{Output}: Parses the motif seeds and their associated 
#'       Transcription Factors (TFs) and writes the result to \code{tfbs_motifs.csv}.
#'   }
#'
#'   The output CSV file contains structured string representations of the 
#'   identified motifs and TFs. The parsed text formats are as follows:
#'   \itemize{
#'     \item \strong{Motif Seeds}: Formatted as 
#'       \code{"{motif.id} {motif.family} [{motif.name}]"}
#'       (e.g., \code{"AsnC; MOTIF 557 AsnC [YwrC - Bacillales] Motif_1"}).
#'     \item \strong{Transcription Factors}: Mapped via TF name and protein model, 
#'       formatted as \code{">{TF.name} {motif.id} Protein:{protein_model}"}
#'       (e.g., \code{">jk0144 1269 Protein:144-jk0144"}).
#'   }
#'
#' @return This function does not return an object to the R environment. 
#'   Its primary purpose is the side effect of generating the 
#'   \code{tfbs_motifs.csv} file in the application working directory.
#'
#' @importFrom seqtoolkit bioseq.fasta bioseq.patterns
#' 
#' @examples
#' # This function is typically executed internally by the GCModeller 
#' # virtual cell framework based on the [@app] decorator, rather than 
#' # called directly by the user.
#' # Example of manual configuration checks:
#' # check_build_module("TRN_network")
[@app "tfbs_motif_scanning"]
const tfbs_motif_scanning = function(app, context) {    
    # run motif site scanning if the TRN_network module build in virtual cell is enabled
    if (check_build_module("TRN_network")) {
        # set the TSS upstream region site fasta file
        let upstream_seq = workfile("make_genbank_proj://upstream_locis.fasta");
        # get the motif database directory that contains multiple meme motif model files
        let motifs_db = file.path(get_config("localdb"), "motifs", get_config("domain"));
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
        #   TF.name motif.id protein_model
        #    ------ ----     --------------
        #   >jk0144 1269 Protein:144-jk0144

        write.csv(motifs, file = outfile, row.names = FALSE);
    }
}

