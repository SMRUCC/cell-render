imports "GenBank" from "seqtoolkit";
imports "bioseq.fasta" from "seqtoolkit";
imports "annotation.genomics" from "seqtoolkit";

#' Extract the genbank source
#' 
#' @param app the current workflow app object
#' @param context the workflow context object
#' 
const extract_gene_table = function(app, context) {
    const src = context$configs$src;
    const gbk = read.genbank(src);
    # get workspace dir path of current workflow app module
    const workdir = WorkflowRender::workspace(app);
    const verbose = as.logical(getOption("verbose"));
    const upstream_size = context$configs$up_len || 150;

    if (verbose) {
        print("app workspace for extract gene table:");
        print(workdir);
    }

    # extract the raw genomics fasta sequence
    const genomics_seq = origin.fasta(gbk);
    const genes = genome.genes(genome = gbk);

    print("get genes table:");
    print([genes]::Synonym);
    print("bp size for parse the gene upstream loci:");
    str(upstream_size);

    const locis = genes 
    |> upstream(length = upstream_size) 
    |> lapply(l -> l, names = [genes]::Synonym)
    |> lapply(function(loci) {
        cut_seq.linear(genomics_seq, loci, doNtAutoReverse = TRUE);
    })
    ;

    # export genomics context elements as feature table.
    write.csv(genes, file = `${workdir}/genes.csv`);
    write.PTT_tabular(gbk, file = `${workdir}/context.txt`);

    # extract sequence data and the gene context data for the
    # downstream transcript regulation network analysis
    write.fasta(genomics_seq, file = `${workdir}/source.fasta`);
    write.fasta(locis, file = `${workdir}/upstream_locis.fasta`);

    # str(app);
    # str(src);
    # str(context);
    invisible(NULL);
}