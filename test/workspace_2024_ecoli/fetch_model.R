require(CellRender);
require(GCModeller);

imports "modeller" from "vcellkit";

let cad_registry = CellRender::open_registry("root", 123456, host = "192.168.3.15");
let template = "G:\biocad_registry\test\Escherichia coli str. K-12 substr. MG1655.txt";
let model = cad_registry |> compile_genbank(gbff = template);

# model 
# |> xml()
# |> writeLines(con = file.path(@dir,"MG1655.xml"))
# ;
write.json_model(model , file =file.path(@dir,"MG1655.json"));