imports "bioseq.fasta" from "seqtoolkit";
imports "bioseq.patterns" from "seqtoolkit";

#' Create the motif TFBS motif data from the gene upstream data
#' 
[@app "tfbs_motif_scanning"]
const tfbs_motif_scanning = function(app, context) {    
    # run motif site scanning if the TRN_network module build in virtual cell is enabled
    if (check_build_module("TRN_network")) {
        # set the TSS upstream region site fasta file
        let upstream_seq = workfile("make_genbank_proj://upstream_locis.fasta");
        # get the motif database directory that contains multiple meme motif model files
        let motifs_db = file.path(get_config("localdb"), "motifs", get_config("domain"));
        # make TFBS site scanning on the TSS upstream region sites
        # search site against the reference motif search.
        let motifs = GCModeller::scan_motifs(
            db = motifs_db, 
            seqs = upstream_seq, 
            workdir = workfile(app, "tfbs_sites"), 
            n_threads = get_config("n_threads"),
            identities_cutoff = 0.9
        );
        let outfile = workfile(app, "tfbs_motifs.csv");

        # seeds 
        # $"{motif.id} {motif.family} [{motif.name}]"
        #
        #         motif.id family [name]
        #              --- ---- -------------------   
        # "AsnC; MOTIF 557 AsnC [YwrC - Bacillales] Motif_1" 

        # TF 
        #
        #   TF.name motif.id protein_model
        #    ------ ----     --------------
        #   >jk0144 1269 Protein:144-jk0144

        write.csv(motifs, file = outfile, row.names = FALSE);
    }
}

