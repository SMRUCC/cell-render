# ==============================================================================
# 1. 环境准备与数据读取
# ==============================================================================

# 加载必要的包
required_packages <- c("Seurat", "ggplot2", "dplyr", "cowplot", "RColorBrewer", "patchwork");

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("F:\\datapool\\16s\\20260223-2000-genomes\\singlecells");

# 读取你的数据文件
# 假设第一行是表头
raw_data <- read.csv("group_phylum.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# ==============================================================================
# 2. 数据预处理 (适配 Seurat 格式)
# ==============================================================================

message("正在处理数据格式...")

# 1) 提取元数据
# 第一列是 assembly id，第二列是 taxonomy (Family)
# 确保第一列名字正确，如果没有名字，暂时命名为 "Assembly_ID"
meta_data <- data.frame(
  row.names = rownames( raw_data), # 行名设为 Assembly ID
  Assembly_ID = rownames( raw_data),
  Family = raw_data[, 1]     # 分组信息
)

# --- 【关键修改】处理缺失值 (NA) ---
# 检查 Family 列是否有 NA，并将 NA 替换为 "Unclassified"
# 这样 DimPlot 就不会报错了
meta_data$Family[is.na(meta_data$Family)] <- "Unclassified"

# 如果有些是空字符串 ""，也建议替换
meta_data$Family[meta_data$Family == ""] <- "Unclassified"

print(head(meta_data));

# 2) 提取纯数值矩阵 (去除前两列)
# 检查是否有非数值列并移除 (保留 EC number 列)
ec_data <- raw_data[, -1];

colnames(ec_data) = sprintf("EC:%s", colnames(ec_data));

# 确保所有列都是数值型，将非数值转换为 NA 后补 0
ec_data <- as.data.frame(lapply(ec_data, function(x) as.numeric(as.character(x))))
ec_data[is.na(ec_data)] <- 0

str(ec_data);

# 3) 构建 Seurat 需要的矩阵
# Seurat 要求：行 = Features (Genes/ECs)，列 = Cells (Genomes)
# 因此我们需要转置
input_matrix <- as.matrix(ec_data)
input_matrix <- t(input_matrix)

# 检查是否有全为 0 的行或列，这会导致 PCA 报错
input_matrix <- input_matrix[rowSums(input_matrix) > 0, ]
# 注意：列通常不会全0，因为基因组至少有一些EC

# ==============================================================================
# 3. 构建 Seurat 对象与分析流程
# ==============================================================================

message("正在构建 Seurat 对象并进行分析...")

# 创建 Seurat 对象
ec_seurat <- CreateSeuratObject(counts = input_matrix, project = "Genome_EC_Embedding", meta.data = meta_data)

# --- 归一化与标准化 ---
# 注意：如果你的嵌入矩阵已经是 TF-IDF 归一化过的，NormalizeData 可以跳过
# 但为了保持数据稳定性，这里执行 LogNormalize
ec_seurat <- NormalizeData(ec_seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# 寻找高变特征 这些是区分不同 Family 的关键代谢酶
ec_seurat <- FindVariableFeatures(ec_seurat, selection.method = "vst", nfeatures = 200, verbose = FALSE)

# 数据缩放
ec_seurat <- ScaleData(ec_seurat, verbose = FALSE)

# PCA 降维
# 这里使用所有 Variable Features 进行 PCA
ec_seurat <- RunPCA(ec_seurat, features = VariableFeatures(object = ec_seurat), verbose = FALSE, npcs = 50)

# UMAP 降维
# dims 决定使用多少个主成分，如果样本量少 (<50)，建议 dims=1:5 或 1:10
# 如果报错 "PCA has not been run"，请检查上面的 RunPCA 是否成功
ec_seurat <- RunUMAP(ec_seurat, dims = 1:10, verbose = FALSE)

# ==============================================================================
# 4. 可视化绘图
# ==============================================================================

print(ec_seurat);

# ------------------------
# 图1: 全局分组 UMAP 图
# 展示不同 Family 在代谢功能空间中的分布
# ------------------------
p_dim <- DimPlot(ec_seurat, reduction = "umap", group.by = "Family", 
                 label = TRUE, repel = TRUE, pt.size = 2, 
                 label.size = 5) + 
  labs(title = "Metabolic Landscape by Family",
       subtitle = "Projected via Seurat UMAP") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right")

ggsave("SC_01_UMAP_Family_Groups.png", p_dim, width = 18, height = 8, dpi = 300)
message("图1 UMAP 分组图已保存")

# ------------------------
# 图2: Feature Plot (特征分布图)
# 展示最具区分度的 Top 4 EC 在 UMAP 上的分布
# ------------------------
top_ecs <- head(VariableFeatures(ec_seurat), 4)

# 绘制 Feature Plot
# min.cutoff 设置为 'q10' 可以去除底部 10% 的噪音，使信号更清晰
p_feature <- FeaturePlot(ec_seurat, features = top_ecs, 
                         min.cutoff = "q10", 
                         cols = c("lightgrey", "darkblue"), 
                         ncol = 2, pt.size = 1.5) + 
  plot_annotation(title = "Top Variable ECs Distribution")

ggsave("SC_02_FeaturePlot_TopECs.png", p_feature, width = 12, height = 10, dpi = 600)
message("图2 Feature Plot 已保存")

# ------------------------
# 图3: Dot Plot (点图)
# 展示各 Family 的特异性代谢指纹
# 点大小 = 携带率，颜色 = 平均丰度
# ------------------------

# 挑选前 10 个高变 EC 进行展示
marker_ecs <- head(VariableFeatures(ec_seurat), 30)

p_dot <- DotPlot(ec_seurat, features = marker_ecs, group.by = "Family") + 
  RotatedAxis() + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Specific EC Markers per Family",
       x = "EC Number", y = "Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        panel.grid.major = element_line(color = "grey90"))

ggsave("SC_03_DotPlot_Markers.png", p_dot, width = 10, height = 18, dpi = 600)
message("图3 Dot Plot 已保存")

# ------------------------
# 图4: 小提琴图
# 直观展示某个 EC 在不同 Family 间的差异
# ------------------------
p_vln <- VlnPlot(ec_seurat, features = marker_ecs[1:2], group.by = "Family", pt.size = 0.5, log = TRUE) + 
  NoLegend() +
  plot_annotation(title = "Expression Distribution of Key ECs")

ggsave("SC_04_ViolinPlot.png", p_vln, width = 30, height = 8, dpi = 300)
message("图4 Violin Plot 已保存")

# ------------------------
# 图5: 热图
# 展示高变 EC 在不同基因组中的表达模式
# ------------------------
p_heat <- DoHeatmap(ec_seurat, features = marker_ecs, group.by = "Family", size = 3, angle = 45) + 
  NoLegend()

ggsave("SC_05_Heatmap.png", p_heat, width = 24, height = 9, dpi = 300)
message("图5 Heatmap 已保存")

# ==============================================================================
# 结束
# ==============================================================================
message("所有分析完成！请查看工作目录下的 PNG 图片。")
