#' Single-Cell Style Analysis of Genome Metabolic Embeddings
#'
#' Performs a single-cell analysis-inspired grouping and visualization of
#' genome metabolic embedding data. The function iterates over three
#' taxonomic ranks (Phylum, Class, Order), assigns each genome to its
#' group at that rank, and delegates to \code{\link{singlecells_viz}} for
#' Seurat-based visualization at each level.
#'
#' This approach treats each genome as a "cell" and each EC number
#' metabolic capability as a "gene", enabling the use of single-cell
#' RNA-seq analysis tools for exploring metabolic diversity patterns
#' across microbial communities.
#'
#' @param embedding_file Character. The file path to the metabolic embedding
#'   CSV file (e.g., produced by \code{\link{diamond_embedding}} or
#'   \code{\link[umap]{analysis}}). The file must contain a \code{taxonomy}
#'   column with standard BIOM-format taxonomy strings.
#' @param workdir Character. The working directory where the grouped analysis
#'   results and visualizations will be saved. Defaults to \code{"./"}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of creating the following directory structure under \code{workdir}:
#'   \itemize{
#'     \item \code{singlecells/group_phylum/} — Phylum-level analysis and plots.
#'     \item \code{singlecells/group_class/} — Class-level analysis and plots.
#'     \item \code{singlecells/group_order/} — Order-level analysis and plots.
#'   }
#'   Each subdirectory contains:
#'   \itemize{
#'     \item \code{group_data.csv} — The embedding data with taxonomy column
#'       replaced by the group label at the corresponding rank.
#'     \item Seurat-based visualization PNG files (see \code{\link{singlecells_viz}}).
#'   }
#'
#' @details
#' The analysis pipeline for each taxonomic rank:
#' \enumerate{
#'   \item Reads the embedding CSV file and parses the taxonomy strings.
#'   \item Replaces the \code{taxonomy} column with the taxonomic name at
#'     the current rank (e.g., for rank = "Phylum", replaces with the
#'     phylum-level name).
#'   \item Creates a rank-specific output directory.
#'   \item Delegates to \code{\link{singlecells_viz}} for Seurat-based
#'     clustering and visualization.
#'   \item Exports the grouped data as a CSV file.
#' }
#'
#' @seealso \code{\link{singlecells_viz}} for the Seurat-based visualization
#'   function, \code{\link{diamond_embedding}} for generating the input
#'   embedding data, \code{\link[umap]{analysis}} for UMAP-based embedding.
#'
#' @examples
#' \dontrun{
#' singlecells_analysis(
#'   embedding_file = "results/embedding/umap.csv",
#'   workdir = "results/singlecells"
#' )
#' }
#'
#' @export
const singlecells_analysis = function(embedding_file, workdir = "./") {
    let data = read.csv(embedding_file, row.names = 1, check.names = FALSE);
    let taxon = biom_string.parse(data$taxonomy);
    let group_dir = NULL;

    for(rank in c("Phylum","Class","Order")) {
        data[, "taxonomy"] <- taxonomy_name(taxon, rank = rank) ;
        group_dir = file.path(workdir, "singlecells", `group_${tolower(rank)}`);

        dir.create(group_dir, showWarnings=FALSE);
        native_r(singlecells_viz, list(
            rawdata = data, 
            outputdir = group_dir
        ));

        write.csv(data, file = file.path(group_dir, "group_data.csv"));
    }
}