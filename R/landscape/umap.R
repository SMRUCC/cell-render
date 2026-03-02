const analysis = function(embedding_file, knn = 200, workdir = "./") {
    let data = read.csv(embedding_file, row.names = 1, check.names = FALSE);
    let taxon = biom_string.parse(data$taxonomy);

    data[,"taxonomy"] = NULL;
    workdir           = workdir || dirname(embedding_file);

    print("view of the scientific names:");
    print(taxonomy_name(taxon, rank = "NA"));

    let result = umap(data, dimension = 3, numberOfNeighbors = knn, localConnectivity = 1, method="Cosine");

    result = as.data.frame(result$umap, labels = result$labels);
    result = cbind(result, data[, c("1.-.-.-","2.-.-.-","3.-.-.-","4.-.-.-","5.-.-.-","6.-.-.-","7.-.-.-")]);

    result[, "scientific_name"] = taxonomy_name(taxon, rank = "NA");
    result[, "kingdom"]         = taxonomy_name(taxon, rank = "Kingdom");
    result[, "phylum"]          = taxonomy_name(taxon, rank = "Phylum");
    result[, "class"]           = taxonomy_name(taxon, rank = "Class");
    result[, "order"]           = taxonomy_name(taxon, rank = "Order");
    result[, "family"]          = taxonomy_name(taxon, rank = "Family");
    result[, "genus"]           = taxonomy_name(taxon, rank = "Genus");
    result[, "species"]         = taxonomy_name(taxon, rank = "Species");    

    write.csv(result, file = file.path(workdir, "umap.csv"));

    native_r(genome_scatter_viz, 
        args = list(data = result, outputdir = workdir)
    );
}