require(GCModeller);

imports "compiler" from "vcellkit";
imports "BioCyc" from "annotationKit";

setwd("F:\\ecoli");

let vcell = "F:\\ecoli\\29.0"
|> open.biocyc()
|> compile_biocyc()
;

vcell 
|> xml()
|> writeLines(con = "ecoli.XML")
;