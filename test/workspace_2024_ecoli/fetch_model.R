require(CellRender);

let cad_registry = CellRender::open_registry("root", 123456, host = "192.168.2.233");
let template = "G:\biocad_registry\test\Escherichia coli str. K-12 substr. MG1655.txt";
let model = cad_registry |> compile_genbank(gbff = template);

model 
|> xml()
|> writeLines(con = file.path(@dir,"MG1655.xml"))
;
