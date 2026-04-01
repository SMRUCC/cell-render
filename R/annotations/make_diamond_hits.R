[@app "make_diamond_hits"]
const make_diamond_hits = function(app, context) {
    let diamond = get_config("diamond");
    let local_db = get_config("localdb");
    let proteins = workfile("make_genbank_proj://proteins.fasta");
    let n_threads = get_config("n_threads");
    let workdir = getwd();
    let diamond_blastp = function(db, output) {
        system2(diamond, c("blastp",
            "--db", db, 
            "--query", proteins, 
            "--out", output, 
            "--threads", n_threads
        ))
        ;
    }

    setwd(WorkflowRender::workspace("make_diamond_hits"));

    make_diamond(local_db, diamond);

    diamond_blastp("ec_number", "ec_number.m8");
    diamond_blastp("subcellular", "subcellular.m8");
    diamond_blastp("transcript_factor","transcript_factor.m8");

    setwd(workdir);
}

#' make the diamond database
const make_diamond = function(local_db, diamond = Sys.which("diamond")) {
    let enzyme_db = file.path(local_db,"ec_numbers.fasta");
    let cc_location = file.path(local_db, "subcellular.fasta");
    let tf_db = file.path(local_db, "TF.fasta");

    diamond <- unlist(diamond);

    system2(diamond, c("makedb","--in",enzyme_db, "--db", "ec_number"));
    system2(diamond, c("makedb","--in",cc_location, "--db", "subcellular"));
    system2(diamond, c("makedb","--in",tf_db, "--db", "transcript_factor"));
}