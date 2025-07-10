require(GCModeller);

imports "rawXML" from "vcellkit";
imports "geneExpression" from "phenotype_kit";

setwd(@dir);

let rawdata = open.vcellPack(file = "./result.vcellPack", mode = "read");
let metabolome = time.frames(rawdata, module="Metabolite", symbol_name= TRUE);

metabolome |> write.expr_matrix(file ="./metabolites.csv",
                                id = "molecules");