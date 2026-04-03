imports "GenBank" from "seqtoolkit";
imports "bioseq.fasta" from "seqtoolkit";
imports "annotation.genomics" from "seqtoolkit";

#' Helper function for extract the annotation features from genbank file
#' 
const extract_gbff = function(src, workdir = "./", 
                              upstream_size = 150, 
                              tag_genbank_accid = FALSE, 
                              verbose = TRUE) {

    # load genbank assembly file from a given file path
    let gbk = read.genbank(src);
    # extract the raw genomics fasta sequence
    let genomics_seq = origin.fasta(gbk);
    # extract the gene features from the genbank assembly object
    let genes = genome.genes(genome = gbk);
    let gene_ids = [genes]::Synonym;
    let genbank_id = gbk |> GenBank::accession_id();

    print(`target genome genbank accession id: ${genbank_id}.`);
    print("get genes table:");
    print(gene_ids);
    print("bp size for parse the gene upstream loci:");
    str(upstream_size);

    let locis = genes 
    # extract gene TSS upstream region and then assign the list name with gene ids
    |> upstream(length = upstream_size, is_relative_offset = TRUE) 
    |> lapply(l => l, names = gene_ids)
    |> tqdm()
    # cast each gene TSS upstream location as site fasta sequence
    # by cut site sequence from genomics sequence via the
    # given TSS upstream location data
    |> lapply(function(loci, i) {
        let fa = cut_seq.linear(genomics_seq, loci, nt_auto_reverse = TRUE);
        let id = gene_ids[i];
        let source_tag = {
            if (tag_genbank_accid) {
                genbank_id;
            } else {
                [fa]::Headers;
            }
        }

        # tag the corresponding gene id to the
        # loci site headers
        fasta.headers(fa) <- append(id, source_tag);
        fa;
    })
    ;

    if (verbose) {
        print("view upstream locis:");
        print(fasta.titles(locis));
    }

    # export genomics context elements as feature table.
    write.csv(genes, file = `${workdir}/genes.csv`);
    # export PTT genomics context tabular file
    write.PTT_tabular(gbk, file = `${workdir}/context.txt`);

    # extract sequence data and the gene context data for the
    # downstream transcript regulation network analysis
    write.fasta(genomics_seq, file = `${workdir}/source.fasta`);
    write.fasta(locis, file = `${workdir}/upstream_locis.fasta`);

    invisible(NULL);
}