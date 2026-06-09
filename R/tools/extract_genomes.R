#' Extract Protein Sequences and Gene Features from Batch GenBank Files
#'
#' @description
#' This function processes multiple GenBank files located in a source directory.
#' It extracts protein sequences and gene feature tables from each file and
#' exports them to a specified output directory. Specifically, it generates
#' FASTA files (\code{.faa}) for protein sequences and CSV files (\code{.csv}) for
#' gene annotation tables.
#'
#' @param src A character string specifying the directory path containing the
#'   GenBank source files. The function scans for files with extensions
#'   \code{.gb}, \code{.gbk}, or \code{.gbff}.
#' @param outputdir A character string specifying the directory path where the
#'   extracted data files will be saved. If the directory does not exist, it
#'   should be created beforehand (or the function may fail depending on
#'   the file system permissions).
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following files for each GenBank file found:
#'   \describe{
#'     \item{\code{<accession_id>.faa}}{Protein FASTA file with locus-tag-only
#'       headers.}
#'     \item{\code{<accession_id>.csv}}{Gene feature annotation table in CSV
#'       format.}
#'   }
#'
#' @details
#' For each GenBank file found in \code{src}:
#' \enumerate{
#'   \item Reads and parses the GenBank file.
#'   \item Extracts the GenBank accession ID.
#'   \item Exports protein sequences with locus-tag-only FASTA headers
#'     (using the template \code{"<locus_tag>"}).
#'   \item Exports the gene feature table as a CSV file with ORF annotations.
#' }
#'
#' @seealso \code{\link{extract_gbff}} for extracting upstream sequences
#'   from a single GenBank file, \code{\link{make_genbank_proj}} for the
#'   full project initialization workflow.
#'
#' @examples
#' \dontrun{
#' extract_genomes(
#'   src = "data/genbank_assemblies/",
#'   outputdir = "results/extracted_proteins/"
#' )
#' }
#'
#' @export
const extract_genomes = function(src, outputdir) {
    # scan all genbank source files inside the given source folder 
    for(let file in tqdm(list.files(src, pattern = c("*.gb","*.gbk","*.gbff")))) {
        # read the genbank file and then make component data exports
        let gb = read.genbank(file);
        let genbank_id = gb |> GenBank::accession_id();
        # export the protein fasta sequence data with fasta header
        # title in format of string template only contains the
        # gene locus tag in the fasta header.
        let prots = gb |> protein_seqs(title = "<locus_tag>");
        let genes = as_tabular(gb, ORF = TRUE);
        
        # export the gene annotation result data from genbank file
        # and also the genomics protein fasta sequence data
        write.csv(genes, file = file.path(outputdir, `${genbank_id}.csv`), silent = TRUE);
        write.fasta(prots, file = file.path(outputdir, `${genbank_id}.faa`));
    }

    invisible(NULL);
}