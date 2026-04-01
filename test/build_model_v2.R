require(CellRender);

setwd(@dir);

modelling_cellgraph(src = "Escherichia coli str. K-12 substr. MG1655.gbff", outputdir = "./test_work", 
                                     up_len = 150, 
                                     localdb = NULL, 
                                     diamond = Sys.which("diamond"), 
                                     n_threads = 32);