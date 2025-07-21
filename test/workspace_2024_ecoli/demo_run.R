require(GCModeller);
require(igraph);

imports ["rawXML", "simulator", "modeller"] from "vcellkit";
imports "analysis" from "vcellkit";
imports "debugger" from "vcellkit";

setwd(@dir);

sink(file = "./run.log");

let modelfile  = "F:\\ecoli\\ecoli.xml";
let model      = as.object(read.vcell(modelfile));
let time.ticks = 300;

print("Run virtual cell model:");
print(model);

let vcell = vcell.model(model,unit.test =FALSE);
let mass  = vcell.mass.index(vcell);
let flux  = vcell.flux.index(vcell);

let dynamics = dynamics.default() |> as.object;

dynamics$transcriptionBaseline   = 2;
dynamics$transcriptionCapacity   = 5;
dynamics$productInhibitionFactor = 0.00000125;

dynamics$translationCapacity = 1;
dynamics$proteinMatureBaseline = 2;
dynamics$proteinMatureCapacity = 10;

print("Using dynamics parameter configuration:");
print(dynamics);

let rawXml = "./result.vcellPack";

let engine = vcell
|> engine.load(	
	inits = mass0(model, unit.test = FALSE, random = [100,5000], map = "metacyc") |> set_status( 
						 Intracellular = list(A = 120, B = 10, C = 0),
						 Extracellular = list(A = 1200,   B = 0,  C = 0)
	),
	iterations       = time.ticks, 
	time_resolutions = 50, 	
	showProgress     = TRUE,
	debug            = FALSE,
	unit.test        = FALSE
) 
|> as.object()
;

debugger::dump_core(engine, file = "./core0.txt");

using xml as open.vcellPack(file  = rawXml, mode  = "write", vcell = engine, graph_debug= FALSE) {
	debugger::set_symbols(xml, model);

	print(rawXml);

	# run virtual cell simulation and then 
	# save the result snapshot data files into 
	# target data directory
	engine$AttachBiologicalStorage(xml);
	engine$Run();

	debugger::dump_core(engine, file = "./core1.txt");
}

sink();