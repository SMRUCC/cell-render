const pfam_diamond = function(proteins, workdir = "./", diamond = Sys.which("diamond")) {
    let pfam = file.path(@datadir, "Pfam-A.fas");
    let ws = getwd();

    workdir = normalizePath(workdir);
    diamond = unlist(diamond);

    dir.create(workdir);
    setwd(workdir);
    system2(diamond, c("makedb","--in",proteins, "--db", "target_proteins"), shell=TRUE);
    system2(diamond, c("blastp",
        "-d","target_proteins.dmnd",
        "-q", pfam, 
        "-o","domains_vs_target.tsv",
        "--ultra-sensitive",
        "--matrix","PAM30",
        "--gapopen","9",
        "--gapextend","1",
        "--evalue","10",
        "--masking","0",
        "--comp-based-stats","0",
        "--outfmt","6","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore",
        "-p","24"));

    setwd(ws);
}