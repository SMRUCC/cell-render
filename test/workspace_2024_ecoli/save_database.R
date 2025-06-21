require(GCModeller);

imports "rawXML" from "vcellkit";
imports "geneExpression" from "phenotype_kit";

let rawdata = open.vcellPack(file = relative_work( "./result.vcellPack"), mode = "read");
let metabolome = time.frames(rawdata, module="Metabolites", symbol_name= FALSE);
let cad_lab = open_cadlab(user="root",passwd=123456, host = "biocad.cloud",port = 3306);

cad_lab |> save_expression(exp_id = "Test-Ecoli-20250621", dynaimics = geneExpression::tr( metabolome));