const sequential_batch = function(src, outputdir = "./", args = list()) {
    let workdir = file.path(outputdir, "tmp");
    let release_dir = get_config("release");    
    let genbank_files = list.files(get_config("src"), 
                                pattern = c("*.gb","*.gbk","*.gbff"), 
                                recursive = TRUE);
    
    message(`Build virtualcell community model based ${length(genbank_files)} genbank source files!`);

    for(let file in genbank_files) {
        file |> make_genbank_proj_file(
            release_dir = release_dir,
            workdir = workdir,
            batch_process = TRUE
        );
    }
}