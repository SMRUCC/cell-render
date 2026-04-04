#' try to get the model genomics accession id if there is plasmid data inside the genbank file
const model_accession_id = function(replicons) {
    let accession_id = NULL;

    for(rep in replicons) {
        if (!is.plasmid(rep)) {
            return(GenBank::accession_id(rep));
        } 

        accession_id <- GenBank::accession_id(rep);
    }

    # all data inside this genbank file is plasmid???
    return(accession_id);
}