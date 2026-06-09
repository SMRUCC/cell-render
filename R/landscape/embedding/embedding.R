require(GCModeller);
require(umap);
require(igraph);

imports ["annotation.workflow", "annotation.terms"] from "seqtoolkit";
imports "taxonomy_kit" from "metagenomics_kit";
imports "OTU_table" from "metagenomics_kit";

#' Generate Metabolic Embedding from DIAMOND BLASTP Results
#'
#' Processes DIAMOND BLASTP alignment output to produce a metabolic
#' capability embedding for a set of genomes. The pipeline extracts
#' EC number annotation terms from BLASTP hits, constructs genome-level
#' metabolic feature vectors, applies TF-IDF vectorization with
#' hierarchical aggregation, and exports the resulting embedding matrix.
#'
#' @param diamond_result Character. The file path to the DIAMOND BLASTP
#'   result file (in \code{.m8} tabular format) or a directory containing
#'   multiple \code{.m8} files.
#' @param workdir Character. The working directory where intermediate and
#'   final output files will be written. Defaults to \code{"./"}.
#' @param union_contigs Integer. The minimum number of genomes in which an
#'   EC number must appear to be included in the final embedding. This
#'   parameter controls the feature selection threshold for the TF-IDF
#'   vectorizer. Defaults to 250.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following files to \code{workdir}:
#'   \itemize{
#'     \item \code{ec_terms.csv} — Extracted EC number annotation terms per genome.
#'     \item \code{genomes_metabolic.jsonl} — Genome metabolic feature vectors
#'       in JSON Lines format.
#'     \item \code{metabolic_embedding.csv} — The final TF-IDF embedding matrix.
#'   }
#'
#' @details
#' The embedding generation pipeline:
#' \enumerate{
#'   \item \strong{Term extraction}: Reads the DIAMOND result stream and
#'     extracts EC number-based metabolic annotation terms via
#'     \code{m8_metabolic_terms()}.
#'   \item \strong{Feature vector construction}: Converts the term stream
#'     into genome-level metabolic feature vectors using
#'     \code{make_vectors()}.
#'   \item \strong{JSONL export}: Writes the feature vectors in JSON Lines
#'     format for interoperability.
#'   \item \strong{TF-IDF vectorization}: Applies TF-IDF transformation
#'     with hierarchical aggregation, filtering EC numbers that appear
#'     in fewer than \code{union_contigs} genomes.
#'   \item \strong{Embedding export}: Writes the final embedding matrix
#'     as a CSV file.
#' }
#'
#' @seealso \code{\link{pangenome_analysis}} for the upstream pan-genome
#'   analysis that produces the DIAMOND results,
#'   \code{\link[umap]{analysis}} for UMAP dimensionality reduction on
#'   the embedding.
#'
#' @examples
#' \dontrun{
#' diamond_embedding(
#'   diamond_result = "results/blastp/ec_number.m8",
#'   workdir = "results/embedding",
#'   union_contigs = 100
#' )
#' }
#'
#' @export
const diamond_embedding = function(diamond_result, workdir = "./", union_contigs = 250) {
    let rawdata = read_m8(diamond_result, stream = TRUE);
    let ec_terms = file.path(workdir, "ec_terms.csv");
    let metabolic_data  = file.path(workdir, "genomes_metabolic.jsonl");
    let ec_emebedding = file.path(workdir, "metabolic_embedding.csv");
    let stream = open.stream(ec_terms, type = "terms");

    stream.flush(m8_metabolic_terms(rawdata), stream);
    stream = open.stream(ec_terms,type = "terms", ioRead = TRUE);

    let genomes = make_vectors(stream, stream = TRUE);

    write_genomes_jsonl(genomes, file = metabolic_data);

    let models = read.jsonl(file = metabolic_data, what = "genome_vector");
    let vec = models |> tfidf_vectorizer(union_contigs = union_contigs, hierarchical  = TRUE);

    write.csv(vec, file = ec_emebedding );
}


