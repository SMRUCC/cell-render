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

    write.fasta(genomics_seq, file = `${workdir}/source.fasta`);
    write.csv(genes, file = `${workdir}/genes.csv`);

    str(app);
    str(src);
    str(context);
}