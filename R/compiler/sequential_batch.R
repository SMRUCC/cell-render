const sequential_batch = function(src, outputdir = "./", args = list()) {
    let workdir = file.path(outputdir, "tmp");
    let release_dir = get_config("release");    
    let genbank_files = list.files(get_config("src"), 
                                pattern = c("*.gb","*.gbk","*.gbff"), 
                                recursive = TRUE)
    ;    
    let enable_blastp_cache = as.logical(get_config("enable_blastp_cache"));
    # reference database dir path
    let localdb = get_config("localdb");  
    let diamond_blastp = diamond_interop(); 

    message(`Build virtualcell community model based ${length(genbank_files)} genbank source files!`);
    setwd(workdir);
    # make reference database
    make_diamond(localdb, diamond);

    for(let file in genbank_files) {
        let model_proj = file |> make_genbank_proj_file(
            release_dir = release_dir,
            workdir = workdir,
            batch_process = TRUE
        );
        let model_id = basename(model_proj );
        let proteins = file.path(model_proj , "proteins.fasta");
        let protein_pfam = file.path(dirname(proteins), "Pfam.csv");

        model_dir <- file.path(workdir, model_id);

        # create workspace dir for save diamond blastp result
        dir.create(model_dir);

        message(`make search for: ${proteins}`);
        message(`diamond blastp export to: ${model_dir}`);

        let ec_out = file.path(model_dir, "ec_number.m8");
        let cc_out = file.path(model_dir, "subcellular.m8");
        let tf_out = file.path(model_dir, "transcript_factor.m8");
        let check_size = 4[KB];
        
        if (file.size(protein_pfam) < check_size) {
            # no cache data
            # run process
            pfam_diamond(
                proteins, 
                workdir = dirname(proteins), 
                diamond = diamond
            );
            
        }        
        
        let check_cache =  (file.size(ec_out) > check_size) 
                        && (file.size(cc_out) > check_size) 
                        && (file.size(tf_out) > check_size)
        ;

        check_cache <- enable_blastp_cache && check_cache;


        if (!check_cache ) {                
            # then run diamond blastp search against the reference database
            diamond_blastp("ec_number", proteins, output = ec_out);
            diamond_blastp("subcellular", proteins, output = cc_out);
            diamond_blastp("transcript_factor",proteins, output = tf_out);
        }

        message(`[${model_id}] diamond blastp search job done!`);

        
    }
}