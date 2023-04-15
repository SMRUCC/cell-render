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
    const workdir = WorkflowRender::workspace(app);

    # extract the raw genomics fasta sequence
    const genomics_seq = origin.fasta(gbk);
    const genes = genome.genes(genome = gbk);

    # extract sequence data and the gene context data for the 
    # downstream transcript regulation network analysis
    write.fasta(genomics_seq, file = `${workdir}/source.fasta`);
    write.csv(genes, file = `${workdir}/genes.csv`);
    write.PTT_tabular(gbk, file = `${workdir}/context.txt`);

    str(app);
    str(src);
    str(context);
}