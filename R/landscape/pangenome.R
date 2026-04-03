imports "pangenome" from "comparative_toolkit";
imports "annotation.workflow" from "seqtoolkit";
imports "annotation.terms" from "seqtoolkit";

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
    # run protein annotation search via diamond blastp search
    batch_diamond(source_dir, blastp_dir, 
        diamond   = normalizePath(unlist(diamond)), 
        n_threads = n_threads);
    
    # make the pan-genome analysis at here
    let links <- new ortho_groups();

    # scan for all diamond blastp result files
    for(let align_file in tqdm(list.files(blastp_dir, pattern = "*.m8"))) {
        # read the diamond blastp alignment file
        align_file 
        |> read_m8()
        # and then set orthology group for make gene family
        |> assign_terms(top_best = TRUE, filter_unknown = TRUE, identities_cut = 30)
        |> set_ortho_group(uf = links)
        ;
    }

    let genome_files = list.files(source_dir, pattern = "*.csv");
    let bin = pangenome::build_context(
        genomes = as.list(genome_files, names = basename(genome_files)) 
                    |> lapply(file => read_genetable(file))
    );
    # run pan-genome analysis
    let result = bin |> pangenome::analysis(orthologSet = links);
    
    # export the analysis result to the result dir
    writeLines(pangenome::report_html(result), con = file.path(result_dir, "analysis_report.html"));
    write.csv(as.data.frame(pangenome::sv_table(result)), file = file.path(result_dir, "SV_result.csv"));
    write.csv(as.data.frame(pangenome::pav_table(result)), file = file.path(result_dir, "PAV_result.csv"));

    invisible(NULL);
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