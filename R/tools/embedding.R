require(GCModeller);
require(umap);
require(igraph);

imports ["annotation.workflow", "annotation.terms"] from "seqtoolkit";
imports "taxonomy_kit" from "metagenomics_kit";
imports "OTU_table" from "metagenomics_kit";

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

const export_tree = function(embedding_file, workdir = "./", node_equals = 0.999) {
    let data = read.OTUtable(embedding_file);
    let graph = OTU_table::makeTreeGraph(
        data |> OTU_table::impute_missing(), 
        equals = node_equals,
        rank_colors = "Class"
    );

    igraph::save.network(graph, file.path(workdir, "metabolic_tree"));
}

const singlecells_analysis = function(embedding_file, workdir = "./") {
    let data = read.csv(embedding_file, row.names = 1, check.names = FALSE);
    let taxon = biom_string.parse(data$taxonomy);
    let group_dir = NULL;

    for(rank in c("Phylum","Class","Order")) {
        data[, "taxonomy"] <- taxonomy_name(taxon, rank = rank) ;
        group_dir = file.path(workdir, "singlecells", `group_${tolower(rank)}`);

        dir.create(group_dir);
        native_r(singlecells_viz, list(
            rawdata = data, 
            outputdir = group_dir
        ));

        write.csv(data, file = file.path(group_dir, "group_data.csv"));
    }
}

const analysis = function(embedding_file, workdir = NULL) {
    let data = read.csv(embedding_file, row.names = 1, check.names = FALSE);
    let taxon = biom_string.parse(data$taxonomy);

    data[,"taxonomy"]=NULL;

    print("view of the scientific names:");
    print(taxonomy_name(taxon, rank = "NA"));

    let result = umap(data, dimension = 9, numberOfNeighbors = 100, localConnectivity = 1, method="Cosine");

    result = as.data.frame(result$umap, labels = result$labels);
    result[, "name"] = taxonomy_name(taxon, rank = "NA") ;
    result[,"phylum"] = taxonomy_name(taxon, rank = "Phylum") ;
    result[,"class"] = taxonomy_name(taxon, rank = "Class") ;
    result[,"order"] = taxonomy_name(taxon, rank = "Order") ;

    workdir = workdir || dirname(embedding_file);

    write.csv(result, file = file.path(workdir, "umap.csv") );

    native_r(genome_scatter_viz, 
        args = list(data = result, outputdir = workdir)
    );
}

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