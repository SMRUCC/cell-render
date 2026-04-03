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
#'   executable. Defaults to the system PATH via \code{Sys.which("diamond")}.
#' @param n_threads An integer specifying the number of CPU threads to use for 
#'   the DIAMOND BLASTP search. Defaults to 32.
#'
#' @details 
#' The pipeline proceeds through the following steps:
#' \itemize{
#'   \item \strong{Extraction}: Extracts protein FASTA (.faa) and gene annotation 
#'     (.csv) files into a \code{source} subdirectory.
#'   \item \strong{Annotation}: Runs DIAMOND BLASTP against an internal EC number 
#'     reference database using \code{\link{batch_diamond}}.
#'   \item \strong{Orthology Assignment}: Reads BLASTP results (.m8), assigns terms 
#'     (keeping only the top best hit with >30\% identity, filtering unknowns), and 
#'     constructs orthology groups.
#'   \item \strong{Pan-genome Construction}: Builds the pangenome context and 
#'     executes the core/accessory analysis using the defined ortholog sets.
#'   \item \strong{Export}: Saves an interactive HTML report, an SV result table 
#'     (CSV), and a PAV result table (CSV) to \code{result_dir}.
#' }
#'
#' @return Returns \code{invisible(NULL)}. All results are written directly to 
#'   the \code{result_dir} as side effects.
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