require(GCModeller);

let vcell = compile_network(
    A +2*B = 3*C+D,
    A = 2*C + 2 * D,
    B = 2 * C
);

vcell = vcell.model(vcell);

let engine = vcell |> engine.load(	
	inits = mass0(vcell, random = [100,5000]) |> set_status( 
						 Intracellular = list(A = 120, B = 10, C = 0),
						 Extracellular = list(A = 1200,   B = 0,  C = 0)
	),
	iterations       = 100, 
	time_resolutions = 1000, 	
	showProgress     = TRUE
) 
|> set_data_driver(data.frame)
|> run()
;

