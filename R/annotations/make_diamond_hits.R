[@app "make_diamond_hits"]
const make_diamond_hits = function(app, context) {
    let diamond = get_config("diamond");
    let local_db = get_config("localdb");
    let proteins = workfile("make_genbank_proj://proteins.fasta");
    let workdir = getwd();

    setwd(WorkflowRender::workspace("make_diamond_hits"));

    make_diamond(local_db, diamond);

    system2(diamond, c("blastp","--db","ec_number", "--query", proteins, "--out", "ec_number.m8", "--threads", n_threads));
    system2(diamond, c("blastp","--db","subcellular", "--query", proteins, "--out", "subcellular.m8", "--threads", n_threads));
    system2(diamond, c("blastp","--db","transcript_factor", "--query", proteins, "--out", "transcript_factor.m8", "--threads", n_threads));

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