require(CellRender);

let cad_registry = CellRender::open_registry("root", 123456, host = "192.168.3.233");

cad_registry |> reaction_pool(repo = "Z:/zzzzz.hds");