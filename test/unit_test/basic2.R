require(GCModeller);

imports ["simulator" "compiler"] from "vcellkit";

let vcell = compile_network(
    A == 2 * B,
    A == B + 3 * C
);
let engine = vcell.model(vcell) |> engine.load(	
	inits = mass0(vcell, random = [100,5000]) |> set_status( 
						Intracellular = list(A = 0, B = 20, C = 0, D = 100)
	),
	iterations       = 30, 
	time_resolutions = 1000, 	
	showProgress     = TRUE
) 
|> attach_memorydataset()
|> run()
;

let result_mass = as.data.frame([engine]::dataStorageDriver, mass = TRUE);
let result_flux = as.data.frame([engine]::dataStorageDriver, mass = FALSE);

print(result_mass);
print(result_flux);

require(ggplot);

result_mass[,"time.axis"] = as.numeric(rownames(result_mass));
result_flux[,"time.axis"] = as.numeric(rownames(result_flux));

bitmap(file = relative_work("basic2_mass.png"), width = 3600, height = 1920) {
    let p = ggplot(result_mass, padding = "padding: 5% 30% 10% 7%;");

    for(name in colnames(result_mass)) {
        if (name != "time.axis") {
            p = p + geom_line(aes(x = "time.axis", y = name), width = 8);
        }        
    }

    p;
}

bitmap(file = relative_work("basic2_flux.png"), width = 3600, height = 1920) {
    let p = ggplot(result_flux, padding = "padding: 5% 30% 10% 7%;");

    for(name in colnames(result_flux)) {
        if (name != "time.axis") {
            p = p + geom_line(aes(x = "time.axis", y = name), width = 8);
        }        
    }

    p;
}