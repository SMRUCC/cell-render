require(GCModeller);
require(igraph);

imports "OTU_table" from "metagenomics_kit";

setwd(@dir);

let data = read.OTUtable("metabolic_embedding.csv");
let graph = OTU_table::makeTreeGraph(
    data |> OTU_table::impute_missing(), 
    equals = 0.999,
    rank_colors = "Class"
);

igraph::save.network(graph, "./metabolic_tree/");