require(GCModeller);

imports ["rawXML", "simulator", "modeller"] from "vcellkit";

setwd(@dir);

let modelfile  = "./ecoli.XML";
let model      = as.object(read.vcell(path = modelfile));
let time.ticks = 1000;

print("Run virtual cell model:");
print(model);

let vcell = vcell.model(model);
let mass  = vcell.mass.index(vcell);
let flux  = vcell.flux.index(vcell);

let dynamics = dynamics.default() :> as.object;

dynamics$transcriptionBaseline   = 200;
dynamics$transcriptionCapacity   = 500;
dynamics$productInhibitionFactor = 0.00000125;

print("Using dynamics parameter configuration:");
print(dynamics);

let rawXml = "result.vcXML";

let engine = vcell
|> engine.load(	
	inits            = mass0(model),
	iterations       = time.ticks, 
	time_resolutions = 1000, 	
	showProgress     = TRUE
) 
|> as.object()
;

using xml as open.vcellXml(file  = rawXml, mode  = "write", vcell = engine) {
	print(rawXml);

	# run virtual cell simulation and then 
	# save the result snapshot data files into 
	# target data directory
	engine$AttachBiologicalStorage(xml);
	engine$Run();
}