require(CellRender);

let cad_registry = CellRender::open_registry("root", 123456, host = "192.168.3.233");

local_protDb(cad_registry, dbfile = "Z:/aaa.fas");