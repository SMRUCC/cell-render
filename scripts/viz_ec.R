# ============================================================================
# EC Number Hierarchical Embedding Visualization
# 可视化基因组间7大类酶的丰度差异
# ============================================================================

# 安装和加载必要的包
required_packages <- c("tidyverse", "ggplot2", "ggrepel", "pheatmap", 
                       "RColorBrewer", "scales", "patchwork", "ggpubr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# 1. 数据读取和预处理
# ============================================================================

#' 读取TF-IDF嵌入结果
#' @param file_path CSV文件路径，包含层级EC嵌入结果
#' @return 包含genome_id和EC特征的data frame
read_hierarchical_embedding <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 检查是否有必要的列
  if (!"assembly_id" %in% colnames(df)) {
    stop("Error: 'assembly_id' column not found in the data!")
  }
  
  cat("Data loaded successfully!\n")
  cat(sprintf("  Genomes: %d\n", nrow(df)))
  cat(sprintf("  Features: %d\n", ncol(df) - 1))
  
  return(df)
}

#' 提取第一层级EC特征（7大类酶）
#' @param df 完整的嵌入数据框
#' @return 只包含第一层级EC的数据框
extract_level1_ec <- function(df) {
  # 第一层级的EC编号模式
  level1_pattern <- "^EC:[1-7]\\.\\.-\\.\\-$"
  
  # 提取匹配的列名
  level1_cols <- grep(level1_pattern, colnames(df), value = TRUE)
  
  if (length(level1_cols) == 0) {
    stop("No Level-1 EC features found! Please check the column names.")
  }
  
  cat(sprintf("\nFound %d Level-1 EC features:\n", length(level1_cols)))
  print(level1_cols)
  
  # 提取相关列
  level1_df <- df %>%
    select(assembly_id, any_of(c("taxonomy", "species")), all_of(level1_cols))
  
  # 创建简化的列名（去掉 EC: 和 .-.-.-）
  new_names <- gsub("EC:([1-7])\\.\\.-\\.\\-", "EC\\1", level1_cols)
  colnames(level1_df)[colnames(level1_df) %in% level1_cols] <- new_names
  
  # 添加酶类标签
  enzyme_labels <- c(
    "EC1" = "Oxidoreductases",
    "EC2" = "Transferases",
    "EC3" = "Hydrolases",
    "EC4" = "Lyases",
    "EC5" = "Isomerases",
    "EC6" = "Ligases",
    "EC7" = "Translocases"
  )
  
  # 检查是否所有7类都存在
  missing_ecs <- setdiff(names(enzyme_labels), colnames(level1_df))
  if (length(missing_ecs) > 0) {
    warning(sprintf("Missing EC classes: %s", paste(missing_ecs, collapse = ", ")))
  }
  
  # 转换为长格式，方便绘图
  level1_long <- level1_df %>%
    pivot_longer(
      cols = starts_with("EC"),
      names_to = "EC_class",
      values_to = "value"
    ) %>%
    mutate(
      EC_number = as.integer(gsub("EC", "", EC_class)),
      Enzyme_name = enzyme_labels[EC_class],
      Enzyme_name = factor(Enzyme_name, levels = enzyme_labels)
    )
  
  return(list(
    wide = level1_df,
    long = level1_long,
    enzyme_labels = enzyme_labels
  ))
}

# ============================================================================
# 2. 雷达图绘制函数
# ============================================================================

#' 绘制单个基因组的雷达图
#' @param data 宽格式数据（一行一个基因组）
#' @param genome_id 要绘制的基因组ID
#' @param title 图表标题
#' @return ggplot对象
plot_single_radar <- function(data, genome_id, title = NULL) {
  # 筛选特定基因组
  genome_data <- data %>%
    filter(assembly_id == genome_id)
  
  if (nrow(genome_data) == 0) {
    stop(sprintf("Genome '%s' not found in data!", genome_id))
  }
  
  # 提取EC值
  ec_cols <- paste0("EC", 1:7)
  values <- as.numeric(genome_data[1, ec_cols])
  
  # 创建雷达图数据
  radar_df <- data.frame(
    Enzyme = c("Oxidoreductases", "Transferases", "Hydrolases", 
               "Lyases", "Isomerases", "Ligases", "Translocases"),
    Value = values
  )
  
  # 添加首行以闭合雷达图
  radar_df <- rbind(radar_df, radar_df[1, ])
  radar_df$Angle <- seq(0, 360, length.out = nrow(radar_df))
  
  # 绘制雷达图
  p <- ggplot(radar_df, aes(x = Angle, y = Value)) +
    geom_polygon(fill = "steelblue", alpha = 0.3, color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 3) +
    coord_polar(start = -pi/7) +
    scale_x_continuous(
      breaks = seq(0, 360, length.out = 8)[-8],
      labels = c("Oxidoreductases", "Transferases", "Hydrolases", 
                 "Lyases", "Isomerases", "Ligases", "Translocases")
    ) +
    labs(
      title = ifelse(is.null(title), sprintf("Genome: %s", genome_id), title),
      x = "",
      y = "TF-IDF Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "gray90"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

#' 绘制多个基因组的雷达图（重叠显示）
#' @param data 宽格式数据
#' @param genome_ids 要比较的基因组ID向量
#' @return ggplot对象
plot_multiple_radar <- function(data, genome_ids) {
  # 筛选基因组
  selected_data <- data %>%
    filter(assembly_id %in% genome_ids)
  
  if (nrow(selected_data) < length(genome_ids)) {
    warning("Some genomes not found in data!")
  }
  
  # 转换为长格式并准备雷达图数据
  ec_cols <- paste0("EC", 1:7)
  
  radar_list <- lapply(1:nrow(selected_data), function(i) {
    genome_id <- selected_data$assembly_id[i]
    values <- as.numeric(selected_data[i, ec_cols])
    
    df <- data.frame(
      Enzyme = c("Oxidoreductases", "Transferases", "Hydrolases", 
                 "Lyases", "Isomerases", "Ligases", "Translocases"),
      Value = values,
      Genome = genome_id
    )
    df <- rbind(df, df[1, ])  # 闭合
    df$Angle <- seq(0, 360, length.out = nrow(df))
    return(df)
  })
  
  radar_df <- do.call(rbind, radar_list)
  
  # 定义颜色
  n_genomes <- length(unique(radar_df$Genome))
  colors <- colorRampPalette(brewer.pal(min(n_genomes, 12), "Set3"))(n_genomes)
  
  # 绘制
  p <- ggplot(radar_df, aes(x = Angle, y = Value, color = Genome, fill = Genome)) +
    geom_polygon(alpha = 0.2, size = 1) +
    coord_polar(start = -pi/7) +
    scale_x_continuous(
      breaks = seq(0, 360, length.out = 8)[-8],
      labels = c("Oxidoreductases", "Transferases", "Hydrolases", 
                 "Lyases", "Isomerases", "Ligases", "Translocases")
    ) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(
      title = "Comparison of Enzyme Class Abundance",
      x = "",
      y = "TF-IDF Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  return(p)
}

# ============================================================================
# 3. 条形图绘制函数
# ============================================================================

#' 绘制堆叠条形图
#' @param data_long 长格式数据
#' @param top_n 显示前N个基因组（按总丰度排序）
#' @return ggplot对象
plot_stacked_barplot <- function(data_long, top_n = 20) {
  # 计算每个基因组的总丰度并排序
  genome_totals <- data_long %>%
    group_by(assembly_id) %>%
    summarise(total = sum(value, na.rm = TRUE)) %>%
    arrange(desc(total)) %>%
    head(top_n)
  
  # 筛选top N基因组
  plot_data <- data_long %>%
    filter(assembly_id %in% genome_totals$assembly_id) %>%
    mutate(assembly_id = factor(assembly_id, levels = genome_totals$assembly_id))
  
  # 绘制
  p <- ggplot(plot_data, aes(x = assembly_id, y = value, fill = Enzyme_name)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = sprintf("Enzyme Class Distribution (Top %d Genomes)", top_n),
      x = "Genome Assembly ID",
      y = "TF-IDF Weighted Abundance",
      fill = "Enzyme Class"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major.x = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  return(p)
}

#' 绘制分组条形图
#' @param data_long 长格式数据
#' @param genome_ids 要显示的基因组ID
#' @return ggplot对象
plot_grouped_barplot <- function(data_long, genome_ids) {
  # 筛选基因组
  plot_data <- data_long %>%
    filter(assembly_id %in% genome_ids) %>%
    mutate(assembly_id = factor(assembly_id, levels = genome_ids))
  
  # 绘制
  p <- ggplot(plot_data, aes(x = Enzyme_name, y = value, fill = assembly_id)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = "Enzyme Class Comparison Across Genomes",
      x = "Enzyme Class",
      y = "TF-IDF Weighted Abundance",
      fill = "Genome"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}

# ============================================================================
# 4. 热图绘制函数
# ============================================================================

#' 绘制热图
#' @param data_wide 宽格式数据
#' @param top_n 显示前N个基因组
#' @param show_rownames 是否显示行名
#' @return pheatmap对象
plot_heatmap <- function(data_wide, top_n = 30, show_rownames = FALSE) {
  # 提取EC列
  ec_cols <- paste0("EC", 1:7)
  
  # 计算总丰度并排序
  data_wide$total <- rowSums(data_wide[, ec_cols], na.rm = TRUE)
  data_sorted <- data_wide %>%
    arrange(desc(total)) %>%
    head(top_n)
  
  # 准备矩阵
  mat <- as.matrix(data_sorted[, ec_cols])
  rownames(mat) <- data_sorted$assembly_id
  colnames(mat) <- c("Oxidoreductases", "Transferases", "Hydrolases", 
                     "Lyases", "Isomerases", "Ligases", "Translocases")
  
  # 标准化（按行）
  mat_scaled <- t(scale(t(mat)))
  
  # 绘制热图
  p <- pheatmap(
    mat_scaled,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = show_rownames,
    show_colnames = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = sprintf("Enzyme Class Abundance Heatmap (Top %d Genomes)", top_n),
    fontsize = 12,
    fontsize_col = 11,
    angle_col = 45,
    border_color = NA
  )
  
  return(p)
}

# ============================================================================
# 5. 主函数：一键生成所有可视化
# ============================================================================

#' 主可视化函数
#' @param file_path TF-IDF结果文件路径
#' @param output_dir 输出目录
#' @param selected_genomes 用于雷达图比较的基因组ID（可选）
visualize_ec_level1 <- function(file_path, 
                                output_dir = "./ec_visualization",
                                selected_genomes = NULL) {
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("\nCreated output directory: %s\n", output_dir))
  }
  
  # 读取数据
  cat("\n=== Reading Data ===\n")
  df <- read_hierarchical_embedding(file_path)
  
  # 提取第一层级EC
  cat("\n=== Extracting Level-1 EC Features ===\n")
  ec_data <- extract_level1_ec(df)
  
  # 如果没有指定基因组，选择总丰度最高的几个
  if (is.null(selected_genomes)) {
    genome_totals <- ec_data$wide %>%
      mutate(total = rowSums(select(., starts_with("EC")), na.rm = TRUE)) %>%
      arrange(desc(total))
    selected_genomes <- head(genome_totals$assembly_id, 5)
    cat(sprintf("\nAuto-selected top 5 genomes for comparison:\n"))
    print(selected_genomes)
  }
  
  # 1. 绘制单个基因组雷达图
  cat("\n=== Plotting Individual Radar Charts ===\n")
  for (genome_id in selected_genomes[1:min(3, length(selected_genomes))]) {
    p <- plot_single_radar(ec_data$wide, genome_id)
    
    output_file <- file.path(output_dir, sprintf("radar_%s.pdf", genome_id))
    ggsave(output_file, p, width = 8, height = 8, units = "in")
    cat(sprintf("  Saved: %s\n", output_file))
  }
  
  # 2. 绘制多基因组雷达图
  cat("\n=== Plotting Multiple Genome Radar Chart ===\n")
  p_radar_multi <- plot_multiple_radar(ec_data$wide, selected_genomes)
  output_file <- file.path(output_dir, "radar_comparison.pdf")
  ggsave(output_file, p_radar_multi, width = 10, height = 8, units = "in")
  cat(sprintf("  Saved: %s\n", output_file))
  
  # 3. 绘制堆叠条形图
  cat("\n=== Plotting Stacked Barplot ===\n")
  p_stacked <- plot_stacked_barplot(ec_data$long, top_n = 20)
  output_file <- file.path(output_dir, "barplot_stacked.pdf")
  ggsave(output_file, p_stacked, width = 12, height = 6, units = "in")
  cat(sprintf("  Saved: %s\n", output_file))
  
  # 4. 绘制分组条形图
  cat("\n=== Plotting Grouped Barplot ===\n")
  p_grouped <- plot_grouped_barplot(ec_data$long, selected_genomes)
  output_file <- file.path(output_dir, "barplot_grouped.pdf")
  ggsave(output_file, p_grouped, width = 10, height = 6, units = "in")
  cat(sprintf("  Saved: %s\n", output_file))
  
  # 5. 绘制热图
  cat("\n=== Plotting Heatmap ===\n")
  p_heatmap <- plot_heatmap(ec_data$wide, top_n = 30, show_rownames = FALSE)
  output_file <- file.path(output_dir, "heatmap.pdf")
  pdf(output_file, width = 8, height = 10)
  print(p_heatmap)
  dev.off()
  cat(sprintf("  Saved: %s\n", output_file))
  
  # 6. 生成汇总报告
  cat("\n=== Generating Summary Report ===\n")
  
  summary_stats <- ec_data$long %>%
    group_by(Enzyme_name) %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Min = min(value, na.rm = TRUE),
      Max = max(value, na.rm = TRUE)
    )
  
  write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  cat("  Saved: summary_statistics.csv\n")
  
  # 保存处理后的数据
  write.csv(ec_data$wide, file.path(output_dir, "level1_ec_data.csv"), row.names = FALSE)
  cat("  Saved: level1_ec_data.csv\n")
  
  # 打印统计信息
  cat("\n=== Summary Statistics ===\n")
  print(summary_stats)
  
  cat(sprintf("\n✓ All visualizations saved to: %s\n", output_dir))
  
  # 返回数据以便后续分析
  return(ec_data)
}

# ============================================================================
# 6. 使用示例
# ============================================================================

# 方法1：使用默认参数（自动选择高丰度基因组）
# ec_data <- visualize_ec_level1("your_tfidf_results.csv")

# 方法2：指定要比较的基因组
# ec_data <- visualize_ec_level1(
#   file_path = "your_tfidf_results.csv",
#   output_dir = "./my_visualization",
#   selected_genomes = c("genome_001", "genome_002", "genome_003")
# )

# 方法3：分步调用，自定义分析
# df <- read_hierarchical_embedding("your_tfidf_results.csv")
# ec_data <- extract_level1_ec(df)
# 
# # 只画雷达图
# p1 <- plot_single_radar(ec_data$wide, "genome_001")
# p2 <- plot_multiple_radar(ec_data$wide, c("genome_001", "genome_002"))
# 
# # 显示图形
# print(p1)
# print(p2)

cat("\n=== R Script for EC Level-1 Visualization ===\n")
cat("Usage:\n")
cat("  ec_data <- visualize_ec_level1('your_tfidf_results.csv')\n\n")
