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