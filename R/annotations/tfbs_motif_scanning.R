imports "bioseq.fasta" from "seqtoolkit";
imports "bioseq.patterns" from "seqtoolkit";

#' Create the motif TFBS motif data from the gene upstream data
#' 
[@app "tfbs_motif_scanning"]
const tfbs_motif_scanning = function(app, context) {    
    let upstream_seq = workfile("make_genbank_proj://upstream_locis.fasta");
    let motifs_db = file.path(get_config("localdb"), "motifs");
    let motifs = GCModeller::scan_motifs(motifs_db, upstream_seq, 
        workdir = workfile(app, "tfbs_sites"), 
        n_threads = get_config("n_threads")
    );
    let outfile = workfile(app, "tfbs_motifs.csv");

    write.csv(motifs, file = outfile, row.names = FALSE);
}