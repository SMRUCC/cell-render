imports "pangenome" from "comparative_toolkit";
imports "annotation.workflow" from "seqtoolkit";
imports "annotation.terms" from "seqtoolkit";

#' Pan-genome Analysis of a Microbial Community
#'
#' This function orchestrates a complete pan-genome analysis pipeline for a
#' set of microbial genomes. It extracts genomic data from GenBank files,
#' performs protein annotation via DIAMOND BLASTP, assigns orthology groups,
#' and calculates both Structural Variation (SV) and Presence/Absence Variation (PAV).
#'
#' @param src A character string specifying the directory path containing the
#'   input NCBI GenBank assembly files.
#' @param result_dir A character string specifying the directory path where
#'   all intermediate and final analysis results will be exported.
#' @param diamond A character string specifying the file path to the DIAMOND
#'   executable. Defaults to \code{Sys.which("diamond")}.
#' @param n_threads An integer specifying the number of CPU threads to use
#'   for DIAMOND BLASTP searches. Defaults to 32.
#' @param identities_cut An integer specifying the percent identity cutoff
#'   for filtering DIAMOND BLASTP hits. Hits below this threshold are
#'   discarded. Defaults to 30.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item Creates an HTML analysis report (\code{analysis_report.html}).
#'     \item Exports a Structural Variation table (\code{SV_result.csv}).
#'     \item Exports a Presence/Absence Variation table (\code{PAV_result.csv}).
#'     \item Writes intermediate protein FASTA and gene annotation CSV files
#'       for each input genome.
#'   }
#'
#' @details
#' The pan-genome analysis pipeline consists of the following steps:
#' \enumerate{
#'   \item \strong{Genome extraction}: Reads each GenBank file and exports
#'     protein sequences (\code{.faa}) and gene annotation tables (\code{.csv}).
#'   \item \strong{DIAMOND BLASTP}: Runs batch DIAMOND BLASTP against the
#'     EC number reference database for all extracted protein sets.
#'   \item \strong{Orthology assignment}: Parses BLASTP results and assigns
#'     orthology group links based on the identity cutoff.
#'   \item \strong{Pan-genome analysis}: Constructs a pan-genome context from
#'     the gene tables and runs the analysis using the orthology links.
#'   \item \strong{Report generation}: Exports the analysis results as an
#'     HTML report and CSV tables.
#' }
#'
#' @seealso \code{\link{extract_genomes}} for genome data extraction,
#'   \code{\link{batch_diamond}} for batch DIAMOND BLASTP execution.
#'
#' @examples
#' \dontrun{
#' pangenome_analysis(
#'   src = "data/genbank_assemblies",
#'   result_dir = "results/pangenome",
#'   diamond = "/usr/bin/diamond",
#'   n_threads = 16,
#'   identities_cut = 30
#' )
#' }
#'
#' @export
const pangenome_analysis = function(src, result_dir, 
                                    diamond = Sys.which("diamond"), 
                                    n_threads = 32, 
                                    skip_blastp = FALSE) {

    let source_dir = file.path(result_dir, "source");
    let blastp_dir = file.path(result_dir, "blastp");

    if (!as.logical(skip_blastp )) {
        # make export of the genomics protein fasta sequence and
        # gene annotation data as files
        extract_genomes(src, outputdir = source_dir);
        # run protein annotation search via diamond blastp search
        batch_diamond(source_dir, blastp_dir, 
            diamond   = normalizePath(unlist(diamond)), 
            n_threads = n_threads);
    }
    
    # make the pan-genome analysis at here
    let links <- list.files(blastp_dir, pattern = "*.m8");
    # scan for all diamond blastp result files
    links <- as.list(links, names = basename(links));
    links <- lapply(tqdm(links), function(align_file) {        
        align_file |> read_m8()
        # read the diamond blastp alignment file
        # and then set orthology group for make gene family
        |> assign_terms(top_best = TRUE, filter_unknown = TRUE, identities_cut = 30)
    });

    let genome_files = list.files(source_dir, pattern = "*.csv");
    let bin = pangenome::build_context(
        genomes = as.list(genome_files, names = basename(genome_files)) 
                    |> tqdm() 
                    |> lapply(file => read_genetable(file))
    );
    # run pan-genome analysis
    let result = bin |> pangenome::analysis(orthologSet = links);
    
    # export the analysis result to the result dir
    writeLines(pangenome::report_html(result), con = file.path(result_dir, "analysis_report.html"));
    write.csv(as.data.frame(pangenome::sv_table(result)), file = file.path(result_dir, "SV_result.csv"));
    write.csv(as.data.frame(pangenome::pav_table(result)), file = file.path(result_dir, "PAV_result.csv"));

    invisible(NULL);
}

