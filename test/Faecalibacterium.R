require(GCModeller);

imports "compiler" from "vcellkit";
imports "BioCyc" from "annotationKit";

setwd("F:\\ecoli");

let vcell = "F:\\ecoli\\gcf_002734145\\29.0"
|> open.biocyc()
|> compile_biocyc()
;

vcell 
|> xml()
|> writeLines(con = "Faecalibacterium.XML")
;