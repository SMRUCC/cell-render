require(GCModeller);
require(umap);

imports "taxonomy_kit" from "metagenomics_kit";

setwd(@dir);

let data = read.csv("metabolic_embedding.csv", row.names = 1, check.names = FALSE);
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

write.csv(result, file = "./umap.csv" );

bitmap(file = "./umap.png") {
    plot(result$dimension_1, result$dimension_2, fill = "white", class = result[,"phylum"],colors = "paper");
}