imports "annotation.workflow" from "seqtoolkit";
imports "proteinKit" from "seqtoolkit";

#' Run Pfam Domain Analysis Using DIAMOND BLASTP
#'
#' Performs protein domain annotation by running DIAMOND BLASTP search
#' against the Pfam-A database. The function creates a DIAMOND database
#' from the input protein sequences, searches them against the Pfam-A
#' reference, and then analyzes the domain architecture of the hits.
#'
#' This function uses a reverse-search strategy: it builds the DIAMOND
#' database from the query proteins and uses the Pfam-A sequences as
#' the search query. This approach is optimized for annotating a small
#' number of genomes against the large Pfam reference.
#'
#' @param proteins Character. The file path to the input protein FASTA
#'   file (\code{.faa}) to be annotated.
#' @param workdir Character. The working directory where DIAMOND
#'   intermediate files will be created. Defaults to \code{"./"}.
#' @param diamond Character. The file path to the DIAMOND executable.
#'   Defaults to \code{Sys.which("diamond")}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following file:
#'   \describe{
#'     \item{\code{Pfam.csv}}{Pfam domain annotation table, written to the
#'       same directory as the input \code{proteins} file. Contains domain
#'       architecture analysis results from \code{proteinKit::analysis_domains()}.}
#'   }
#'
#' @details
#' The analysis pipeline:
#' \enumerate{
#'   \item Changes the working directory to \code{workdir}.
#'   \item Creates a DIAMOND database from the input protein sequences.
#'   \item Runs DIAMOND BLASTP with the Pfam-A reference as query against
#'     the protein database, using ultra-sensitive mode with PAM30 matrix
#'     and relaxed e-value threshold (10) for domain detection.
#'   \item Parses the BLASTP output in m8 format.
#'   \item Analyzes domain architecture using \code{proteinKit::analysis_domains()}.
#'   \item Restores the original working directory.
#'   \item Exports the domain annotation results to \code{Pfam.csv}.
#' }
#'
#' The DIAMOND search parameters are specifically tuned for Pfam domain
#' detection: PAM30 substitution matrix, gap open penalty of 9, gap
#' extension penalty of 1, and composition-based statistics disabled.
#'
#' @seealso \code{\link{make_diamond_hits}} which calls this function
#'   as part of the annotation workflow, \code{\link{make_blastp_term}}
#'   for parsing other BLASTP results.
#'
#' @examples
#' \dontrun{
#' pfam_diamond(
#'   proteins = "results/proteins.fasta",
#'   workdir = "results/pfam_search",
#'   diamond = "/usr/bin/diamond"
#' )
#' }
#'
#' @keywords internal
#' @export
const pfam_diamond = function(proteins, workdir = "./", diamond = Sys.which("diamond")) {
    let pfam = file.path(@datadir, "Pfam-A.fas");
    let ws = getwd();
    let protein_id = basename(proteins);
    let m8_file = `${protein_id}.tsv`;

    workdir  = normalizePath(workdir);
    proteins = normalizePath(proteins);
    diamond  = unlist(diamond);
    
    dir.create(workdir, showWarnings=FALSE);
    setwd(workdir);
    system2(diamond, c("makedb","--in",proteins, "--db", protein_id), shell=TRUE);
    system2(diamond, c("blastp",
        "-d",`${protein_id}.dmnd`,
        "-q", pfam, 
        "-o", m8_file,
        "-p","24",
        "--ultra-sensitive",
        "--matrix","PAM30",
        "--gapopen","9",
        "--gapextend","1",
        "--evalue","10",
        "--masking","0",
        "--comp-based-stats","0",
        "--outfmt","6","qtitle","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), shell=TRUE);

    pfam = read_m8(m8_file);
    pfam = proteinKit::analysis_domains(pfam);

    setwd(ws);

    write.csv(pfam, file = file.path(dirname(proteins), "Pfam.csv"));
}