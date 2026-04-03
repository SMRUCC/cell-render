#' Extract Protein Sequences and Gene Features from Batch GenBank Files
#'
#' @description
#' This function processes multiple GenBank files located in a source directory.
#' It extracts protein sequences and gene feature tables from each file and
#' exports them to a specified output directory. Specifically, it generates
#' FASTA files (`.faa`) for protein sequences and CSV files (`.csv`) for
#' gene annotation tables.
#'
#' @param src A character string specifying the directory path containing the
#'   GenBank source files. The function scans for files with extensions
#'   `.gb`, `.gbk`, or `.gbff`.
#' @param outputdir A character string specifying the directory path where the
#'   extracted data files will be saved. If the directory does not exist, it
#'   should be created beforehand (or the function may fail depending on
#'   implementation).
#'
#' @return Invisible `NULL`. This function is called for its side effects
#'   (writing files to disk) and does not return a value to the console.
#'
#' @details
#' The function iterates through all matching GenBank files in the source
#' directory. For each file:
#' \itemize{
#'   \item Reads the GenBank record.
#'   \item Extracts the accession ID to name the output files.
#'   \item Exports protein sequences in FASTA format, using the gene locus tag
#'         as the sequence header.
#'   \item Exports gene features (ORFs) as a tabular CSV file.
#' }
#'
#' @seealso
#' [read.genbank], [GenBank::accession_id], [write.fasta]
#'
#' @examples
#' \dontrun{
#' # Define source and output directories
#' source_path <- "path/to/genbank_files"
#' output_path <- "path/to/output_data"
#'
#' # Run the extraction
#' extract_genomes(src = source_path, outputdir = output_path)
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