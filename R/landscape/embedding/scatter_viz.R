#' Generate Scatter Plot Visualization of Genome Embedding Results
#'
#' Creates a publication-quality scatter plot from genome embedding data,
#' coloring points by taxonomic phylum and sizing points inversely proportional
#' to the number of genomes in each phylum group. The plot is saved in both
#' PNG and SVG formats.
#'
#' This function uses native R with \pkg{ggplot2} and \pkg{dplyr} for
#' data processing and visualization.
#'
#' @param data A data frame containing the embedding coordinates with at least
#'   the following columns:
#'   \describe{
#'     \item{dimension_1}{Numeric. The first embedding dimension (x-axis).}
#'     \item{dimension_2}{Numeric. The second embedding dimension (y-axis).}
#'     \item{phylum}{Character. The taxonomic phylum assignment for each genome.}
#'   }
#' @param outputdir Character. The directory path where the output plot files
#'   will be saved. Defaults to \code{"./"}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following files to \code{outputdir}:
#'   \itemize{
#'     \item \code{scatter_plot.png} — High-resolution PNG scatter plot
#'       (600 DPI, 12 x 8 inches).
#'     \item \code{scatter_plot.svg} — SVG vector scatter plot
#'       (600 DPI, 12 x 8 inches).
#'   }
#'
#' @details
#' The visualization applies the following design choices:
#' \itemize{
#'   \item Points are colored by \code{phylum} using the default ggplot2
#'     discrete color palette.
#'   \item Point size is mapped to the subset size of each phylum group,
#'     with a \strong{reversed} scale so that larger groups (more genomes)
#'     are rendered with \emph{smaller} points, reducing visual clutter
#'     from dominant phyla.
#'   \item A global alpha of 0.3 is applied for transparency to handle
#'     overlapping points.
#'   \item The minimal theme is used for a clean, publication-ready appearance.
#' }
#'
#' @seealso \code{\link[umap]{analysis}} which generates the embedding data
#'   consumed by this function, \code{\link{singlecells_viz}} for an
#'   alternative Seurat-based visualization approach.
#'
#' @examples
#' \dontrun{
#' # Assuming 'embedding_data' has columns: dimension_1, dimension_2, phylum
#' genome_scatter_viz(
#'   data = embedding_data,
#'   outputdir = "results/figures"
#' )
#' }
#'
#' @export
let genome_scatter_viz = function(data, outputdir = "./") {
  require(ggplot2);
  require(dplyr);

  # 数据预处理：计算每个 Phylum 分组的大小
  # 我们需要给每一行数据都标记上它所属的 Phylum 组有多少个样本
  let data_processed = data 
  |> group_by(phylum) 
  # 计算每组的数量，并存为新列 subset_size
  |> mutate(subset_size = n())  
  |> ungroup()
  ;

  # 绘制散点图
  let p <- ggplot(data_processed, aes(x = dimension_1, y = dimension_2))
    + geom_point(
      aes(
        color = phylum,      # 颜色映射到 phylum
        size = subset_size   # 大小映射到分组子集大小
      ),
      alpha = 0.3            # 设置全局半透明度为 0.3
    )    
    # 设置大小的映射变换：集合越大，点越小
    # trans = "reverse" 会反转坐标轴，即数据值越大(集合越大)，图形尺寸越小
    + scale_size_continuous(
      name = "Subset Size",  # 图例标题
      trans = "reverse",     # 关键：反转大小逻辑
      range = c(2, 6)        # 设置点的大小范围（最小点，最大点)
    )
    # 添加图表标签和主题
    + labs(
      title = "Scatter Plot Grouped by Phylum",
      x = "Dimension 1",
      y = "Dimension 2",
      # 颜色图例标题
      color = "Phylum"            
    ) 
    + theme_minimal()
  ;

  ggsave(file.path( "scatter_plot.png"), plot = p, 
    width = 12, height = 8, 
    dpi = 600);
  ggsave(file.path("scatter_plot.svg"), plot = p, 
    width = 12, height = 8, 
    dpi = 600);
}