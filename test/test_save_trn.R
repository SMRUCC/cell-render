require(CellRender);

modelling_cellgraph(src = "G:\GMNDesigner-DEMO\Acetoin\genbank", outputdir = "K:\models", 
                                     name = NULL,
                                     up_len = 150, 
                                     localdb = NULL, 
                                     diamond = ("diamond"), 
                                     domain = c("bacteria"),
                                     builds = c("TRN_network","Metabolic_network"),
                                     enable_blastp_cache = FALSE,
                                     enzyme_fuzzy = FALSE,
                                     gems_library_mode = TRUE,
                                     batch_mode = c("batch"),
                                     n_threads = 32, 
                                     debug = c("make_TRN"));