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