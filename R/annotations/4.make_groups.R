require(GCModeller);

imports "taxonomy_kit" from "metagenomics_kit";

setwd(@dir);

let data = read.csv("metabolic_embedding.csv", row.names = 1, check.names = FALSE);
let taxon = biom_string.parse(data$taxonomy);

for(rank in c("Phylum","Class","Order")) {
    data[, "taxonomy"] <- taxonomy_name(taxon, rank = rank) ;

    write.csv(data, file = file.path("singlecells", `group_${tolower(rank)}.csv`));
}