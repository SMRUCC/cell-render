require(CellRender);

const genbank_src = "K:\Xanthomonas_campestris_8004_uid15\genbank\CP000050.1.txt"
const outputdir = "K:\Xanthomonas_campestris_8004_uid15\genbank\CP000050.1_test"

CellRender.modelling_cellgraph(genbank_src, outputdir = outputdir);