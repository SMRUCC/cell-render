imports "taxonomy_kit" from "metagenomics_kit";
imports "bifrost" from "seqtoolkit";

#' Compile a Metagenome-Assembled Genome (MAG) as a Virtual Cell Model
#'
#' Processes a metagenome-assembled genome (MAG) contigs assembly by
#' performing gene prediction with Prodigal, then exporting the predicted
#' genes and proteins in multiple standard bioinformatics formats for
#' downstream analysis and virtual cell model construction.
#'
#' @param contigs A character string or object representing the MAG contigs
#'   assembly data. This is typically a FASTA file path or a parsed contigs
#'   object that can be processed by \code{bifrost::prodigal()}.
#' @param tax_str A character string specifying the taxonomy of this MAG
#'   contigs assembly in GTDB-style format, e.g.:
#'   \code{"d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina sp001714685"}.
#'   If a character string is provided, it will be parsed into a structured
#'   taxonomy object via \code{biom_string.parse()}.
#' @param outputdir Character. The directory path where the output files
#'   will be saved. Defaults to \code{"./"}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following files to \code{outputdir}:
#'   \describe{
#'     \item{\code{gene_predicts.csv}}{Gene prediction results in CSV format.}
#'     \item{\code{gene_predicts.gff3}}{Gene predictions in GFF3 format.}
#'     \item{\code{gene_predicts.fna}}{Nucleotide FASTA of predicted gene sequences.}
#'     \item{\code{protein_predicts.faa}}{Protein FASTA of predicted protein sequences.}
#'   }
#'
#' @details
#' The processing pipeline:
#' \enumerate{
#'   \item Runs \code{bifrost::prodigal()} on the contigs to predict
#'     open reading frames (ORFs) with a minimum ORF length of 90 bp.
#'   \item If \code{tax_str} is a character string, parses it into a
#'     structured taxonomy object.
#'   \item Exports the gene prediction results in four standard formats:
#'     CSV (tabular summary), GFF3 (standard annotation format), FASTA
#'     nucleotide (gene sequences), and FASTA protein (translated sequences).
#' }
#'
#' @seealso \code{\link{make_genbank_proj}} for processing complete GenBank
#'   assemblies, \code{\link{build_project}} for compiling annotations into
#'   virtual cell models.
#'
#' @examples
#' \dontrun{
#' MAGs_model(
#'   contigs = "data/MAGs/bin_001.fasta",
#'   tax_str = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis",
#'   outputdir = "results/MAGs/bin_001"
#' )
#' }
#'
#' @export
const MAGs_model = function(contigs, tax_str, outputdir = "./") {
    # make gene predictions from the MAGs contigs
    let genes = bifrost::prodigal(contigs, min.ORF.len = 90);

    if (is.character(tax_str)) {
        tax_str <- biom_string.parse(tax_str);
    }

    # export result to files
    write.csv(as.data.frame(result), file = file.path(outputdir, "gene_predicts.csv"));
    # export the gene prediction result to GFF3 format
    write.gff3(as.gff3(result), file = file.path(outputdir,"gene_predicts.gff3"));
    # export gene/protein fasta sequence to file
    write.fasta(as.genes(result), file = file.path(outputdir,"gene_predicts.fna"));
    write.fasta(as.proteins(result), file = file.path(outputdir,"protein_predicts.faa"));
}