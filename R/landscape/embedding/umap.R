#' UMAP Dimensionality Reduction and Taxonomy-Enriched Visualization
#'
#' Performs UMAP (Uniform Manifold Approximation and Projection) dimensionality
#' reduction on metabolic embedding data, appends full taxonomic annotations
#' to the result, and generates a scatter plot visualization via
#' \code{\link{genome_scatter_viz}}.
#'
#' The function uses Cosine distance as the similarity metric and produces
#' a 3-dimensional UMAP embedding. The first seven EC number hierarchy levels
#' (EC1 through EC7) are preserved as additional columns in the output for
#' downstream analysis.
#'
#' @param embedding_file Character. The file path to the metabolic embedding
#'   CSV file. The file must contain a \code{taxonomy} column with standard
#'   BIOM-format taxonomy strings, and EC number columns including at least
#'   \code{1.-.-.-} through \code{7.-.-.-} for the top-level EC hierarchy.
#' @param knn Integer. The number of nearest neighbors to use for UMAP
#'   graph construction. Larger values produce more global structure,
#'   while smaller values preserve local structure. Defaults to 200.
#' @param workdir Character. The directory path where output files will be
#'   saved. Defaults to \code{"./"}. If not provided, the directory of
#'   \code{embedding_file} is used.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following files to \code{workdir}:
#'   \itemize{
#'     \item \code{umap.csv} — The UMAP embedding result with appended
#'       taxonomic annotations and EC hierarchy columns.
#'     \item \code{scatter_plot.png} — Scatter plot of the UMAP embedding
#'       (via \code{\link{genome_scatter_viz}}).
#'     \item \code{scatter_plot.svg} — SVG version of the scatter plot.
#'   }
#'
#' @details
#' The analysis pipeline:
#' \enumerate{
#'   \item Reads the embedding CSV and parses the taxonomy strings.
#'   \item Removes the \code{taxonomy} column before running UMAP
#'     (to avoid passing non-numeric data to the algorithm).
#'   \item Runs UMAP with 3 dimensions, Cosine distance metric, and the
#'     specified \code{knn} nearest neighbors.
#'   \item Appends the following taxonomic columns to the result:
#'     \code{scientific_name}, \code{kingdom}, \code{phylum}, \code{class},
#'     \code{order}, \code{family}, \code{genus}, \code{species}.
#'   \item Preserves the top-level EC hierarchy columns
#'     (\code{1.-.-.-} through \code{7.-.-.-}).
#'   \item Exports the enriched result as CSV.
#'   \item Generates a scatter plot via \code{\link{genome_scatter_viz}}.
#' }
#'
#' @seealso \code{\link{diamond_embedding}} for generating the input
#'   embedding data, \code{\link{genome_scatter_viz}} for the scatter
#'   plot visualization, \code{\link{singlecells_analysis}} for an
#'   alternative Seurat-based analysis.
#'
#' @examples
#' \dontrun{
#' analysis(
#'   embedding_file = "results/embedding/metabolic_embedding.csv",
#'   knn = 150,
#'   workdir = "results/umap"
#' )
#' }
#'
#' @export
const analysis = function(embedding_file, knn = 200, workdir = "./") {
    let data = read.csv(embedding_file, row.names = 1, check.names = FALSE);
    let taxon = biom_string.parse(data$taxonomy);

    data[,"taxonomy"] = NULL;
    workdir           = workdir || dirname(embedding_file);

    print("view of the scientific names:");
    print(taxonomy_name(taxon, rank = "NA"));

    let result = umap(data, dimension = 3, numberOfNeighbors = knn, localConnectivity = 1, method="Cosine");

    result = as.data.frame(result$umap, labels = result$labels);
    result = cbind(result, data[, c("1.-.-.-","2.-.-.-","3.-.-.-","4.-.-.-","5.-.-.-","6.-.-.-","7.-.-.-")]);

    result[, "scientific_name"] = taxonomy_name(taxon, rank = "NA");
    result[, "kingdom"]         = taxonomy_name(taxon, rank = "Kingdom");
    result[, "phylum"]          = taxonomy_name(taxon, rank = "Phylum");
    result[, "class"]           = taxonomy_name(taxon, rank = "Class");
    result[, "order"]           = taxonomy_name(taxon, rank = "Order");
    result[, "family"]          = taxonomy_name(taxon, rank = "Family");
    result[, "genus"]           = taxonomy_name(taxon, rank = "Genus");
    result[, "species"]         = taxonomy_name(taxon, rank = "Species");    

    write.csv(result, file = file.path(workdir, "umap.csv"));

    native_r(genome_scatter_viz, 
        args = list(data = result, outputdir = workdir)
    );
}