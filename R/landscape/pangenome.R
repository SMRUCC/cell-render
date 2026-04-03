#' pan-genome analysis of the microbial community
#' 
#' @param src a character vector of the directory path that contains the ncbi genbank assembly files.
#' @param result_dir a character vector of the directory path for export the analysis result files. 
#' @param diamond diamond software file path, for make blastp search
#' 
const pangenome = function(src, result_dir, diamond = Sys.which("diamond"), n_threads = 32) {
    let source_dir = file.path(result_dir, "source");
    let blastp_dir = file.path(result_dir, "blastp");
    
    # make export of the genomics protein fasta sequence and
    # gene annotation data as files
    extract_genomes(src, outputdir = source_dir);
    batch_diamond(source_dir, blastp_dir, 
        diamond   = unlist(diamond), 
        n_threads = n_threads);
    


}

#' make diamond blastp search against the reference enzyme database in batch
#' 
const batch_diamond = function(source_dir, result_dir, 
                               diamond = Sys.which("diamond"), 
                               n_threads = 32) {
    let current_dir = getwd();
    # use the ec_numbers.fasta as reference db
    # which is inside the data dir of current R package
    let local_db = file.path(@datadir, "ec_numbers.fasta");

    setwd(blastp_dir);        
    system2(diamond, c("makedb","--in", local_db, "--db", "ec_number"));

    for(let faa in list.files(source_dir, pattern = "*.faa")) {
        system2(diamond, c("blastp",
            "--db", "ec_number", 
            "--query", faa, 
            "--out", file.path(result_dir, `${basename(faa)}.m8`), 
            "--threads", n_threads
        ))
        ;
    }

    setwd(current_dir);
}