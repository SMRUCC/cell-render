#' Export Metabolic Tree Network Visualization
#'
#' Constructs a tree graph from metabolic embedding data and exports it
#' as a network visualization file. The tree graph is built from an
#' OTU-style table using imputation for missing values and a configurable
#' node equality threshold for branch merging.
#'
#' @param embedding_file Character. The file path to the metabolic embedding
#'   data file (OTU table format) containing the feature matrix for tree
#'   construction.
#' @param workdir Character. The directory path where the exported network
#'   file will be saved. Defaults to \code{"./"}.
#' @param node_equals Numeric. The similarity threshold for merging tree
#'   nodes. Nodes with similarity above this value will be merged into
#'   a single branch. Values closer to 1.0 produce more granular trees,
#'   while lower values produce more collapsed trees. Defaults to 0.999.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effect
#'   of writing a network file named \code{metabolic_tree} to \code{workdir}.
#'   The file format is determined by \code{igraph::save.network()}.
#'
#' @details
#' The tree construction pipeline:
#' \enumerate{
#'   \item Reads the embedding data as an OTU table.
#'   \item Imputes missing values in the OTU table.
#'   \item Constructs a tree graph using \code{OTU_table::makeTreeGraph()},
#'     with nodes colored by the \code{Class} taxonomic rank.
#'   \item Saves the resulting network graph to disk.
#' }
#'
#' @seealso \code{\link{diamond_embedding}} for generating the input
#'   embedding data, \code{\link[umap]{analysis}} for UMAP-based
#'   alternative visualization.
#'
#' @examples
#' \dontrun{
#' export_tree(
#'   embedding_file = "results/embedding/metabolic_embedding.csv",
#'   workdir = "results/tree_viz",
#'   node_equals = 0.995
#' )
#' }
#'
#' @export
const export_tree = function(embedding_file, workdir = "./", node_equals = 0.999) {
    let data = read.OTUtable(embedding_file);
    let graph = OTU_table::makeTreeGraph(
        data |> OTU_table::impute_missing(), 
        equals = node_equals,
        rank_colors = "Class"
    );

    igraph::save.network(graph, file.path(workdir, "metabolic_tree"));
}


