imports "bioseq.fasta" from "seqtoolkit";
imports "bioseq.patterns" from "seqtoolkit";

#' Create the motif TFBS motif data from the gene upstream data
#' 
const tfbs_motif_scanning = function(app, context) {
    const src_dir = WorkflowRender::workspace(extract_gene_table);
    const upstream_seq = read.fasta(`${src_dir}/upstream_locis.fasta`);
    const motifs = find_motifs(upstream_seq);
    const workdir = WorkflowRender::workspace(app);

    motifs 
    |> xml
    |> writeLines(con = `${workdir}/tfbs_motifs.xml`)
    ;
}