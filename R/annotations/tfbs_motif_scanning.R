imports "bioseq.fasta" from "seqtoolkit";
imports "bioseq.patterns" from "seqtoolkit";

#' Create the motif TFBS motif data from the gene upstream data
#' 
[@app "tfbs_motif_scanning"]
const tfbs_motif_scanning = function(app, context) {    
    const upstream_seq = read.fasta(workfile("make_genbank_proj://upstream_locis.fasta"));
    const motifs = find_motifs(upstream_seq);
    const outfile = workfile(app, "tfbs_motifs.xml");

    motifs 
    |> xml
    |> writeLines(con = outfile)
    ;
}