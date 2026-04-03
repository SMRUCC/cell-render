#' extract of the batch list of the proteins sequence and gene features inside genbank data
#' 
#' @param src a directory path that should contains multiple genbank files 
#' @param outputdir a directory path for export the data files from the genbank source files
#'  
const extract_genomes = function(src, outputdir) {
    # scan all genbank source files inside the given source folder 
    for(let file in tqdm(list.files(src, pattern = c("*.gb","*.gbk","*.gbff")))) {
        # read the genbank file and then make component data exports
        let gb = read.genbank(file);
        let genbank_id = gb |> GenBank::accession_id();
        # export the protein fasta sequence data with fasta header
        # title in format of string template only contains the
        # gene locus tag in the fasta header.
        let prots = gb |> protein_seqs(title = "<locus_tag>");
        let genes = as_tabular(gb, ORF = TRUE);
        
        # export the gene annotation result data from genbank file
        # and also the genomics protein fasta sequence data
        write.csv(genes, file = file.path(outputdir, `${genbank_id}.csv`));
        write.fasta(prots, file = file.path(outputdir, `${genbank_id}.faa`));
    }

    invisible(NULL);
}