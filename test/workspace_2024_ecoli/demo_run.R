require(GCModeller);
require(igraph);

imports ["rawXML", "simulator", "modeller"] from "vcellkit";
imports "analysis" from "vcellkit";
imports "debugger" from "vcellkit";

setwd(@dir);

let modelfile  = "./MG1655.xml";
let model      = as.object(read.vcell(path = modelfile));
let time.ticks = 1000;

print("Run virtual cell model:");
print(model);

let vcell = vcell.model(model);
let mass  = vcell.mass.index(vcell);
let flux  = vcell.flux.index(vcell);

let dynamics = dynamics.default() |> as.object;

dynamics$transcriptionBaseline   = 200;
dynamics$transcriptionCapacity   = 500;
dynamics$productInhibitionFactor = 0.00000125;

print("Using dynamics parameter configuration:");
print(dynamics);

let rawXml = "./result.vcellPack";

let engine = vcell
|> engine.load(	
	inits            = mass0(model),
	iterations       = time.ticks, 
	time_resolutions = 1000, 	
	showProgress     = TRUE
) 
|> as.object()
;

debugger::dump_core(engine, file = "./core0.txt");

using xml as open.vcellPack(file  = rawXml, mode  = "write", vcell = engine, graph_debug= FALSE) {
	print(rawXml);

	# run virtual cell simulation and then 
	# save the result snapshot data files into 
	# target data directory
	engine$AttachBiologicalStorage(xml);
	engine$Run();

	debugger::dump_core(engine, file = "./core1.txt");
}