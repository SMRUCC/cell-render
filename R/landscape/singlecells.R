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