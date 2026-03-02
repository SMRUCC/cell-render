const export_tree = function(embedding_file, workdir = "./", node_equals = 0.999) {
    let data = read.OTUtable(embedding_file);
    let graph = OTU_table::makeTreeGraph(
        data |> OTU_table::impute_missing(), 
        equals = node_equals,
        rank_colors = "Class"
    );

    igraph::save.network(graph, file.path(workdir, "metabolic_tree"));
}


