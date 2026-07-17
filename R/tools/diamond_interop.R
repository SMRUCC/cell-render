const diamond_interop = function() {
    let n_threads = get_config("n_threads");
    let diamond   = get_config("diamond");  # diamond path
    # a helper wrapper function of the diamond blastp search 
    # commandline calls
    let diamond_blastp = function(db, proteins, output) {
        system2(diamond, c("blastp",
            "--db", db, 
            "--query", proteins, 
            "--out", output, 
            "--threads", n_threads,
            "--outfmt","6", 
            "qseqid","stitle","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"
        ), shell=TRUE)
        ;
    }

    return(diamond_blastp );
}