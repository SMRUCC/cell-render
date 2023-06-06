require(GCModeller);

imports "compiler" from "vcellkit";
imports "BioCyc" from "annotationKit";
imports "GenBank" from "seqtoolkit";

let genome_src = "E:\GCA_000005845.2_ASM584v2_genomic\source\NC_000913.gb";
let vcell = "E:\UniProt\BioCyc\tier1\ecoli\25.1"
|> open.biocyc()
|> compile.biocyc(genomes = read.genbank(genome_src, repliconTable = TRUE))
;

vcell 
|> xml()
|> writeLines(con = `${@dir}/ecoli.XML`)
;