imports "annotation.workflow" from "seqtoolkit";
imports "proteinKit" from "seqtoolkit";

const pfam_diamond = function(proteins, workdir = "./", diamond = Sys.which("diamond")) {
    let pfam = file.path(@datadir, "Pfam-A.fas");
    let ws = getwd();
    let protein_id = basename(proteins);
    let m8_file = `${protein_id}.tsv`;

    workdir  = normalizePath(workdir);
    proteins = normalizePath(proteins);
    diamond  = unlist(diamond);
    
    dir.create(workdir);
    setwd(workdir);
    system2(diamond, c("makedb","--in",proteins, "--db", protein_id), shell=TRUE);
    system2(diamond, c("blastp",
        "-d",`${protein_id}.dmnd`,
        "-q", pfam, 
        "-o", m8_file,
        "-p","24",
        "--ultra-sensitive",
        "--matrix","PAM30",
        "--gapopen","9",
        "--gapextend","1",
        "--evalue","10",
        "--masking","0",
        "--comp-based-stats","0",
        "--outfmt","6","qtitle","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"), shell=TRUE);

    pfam = read_m8(m8_file);
    pfam = proteinKit::analysis_domains(pfam);

    setwd(ws);

    write.csv(pfam, file = file.path(dirname(proteins), "Pfam.csv"));
}