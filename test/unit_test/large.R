require(GCModeller);

imports ["simulator" "compiler"] from "vcellkit";

let vcell = compile_network(
    D_galactose + O2 == D_galacto_hexodialdose + H2O2,
    D_galactono_1_5_lactone + H2O2 == D_galactose + O2,
    D_Galactono_1_5_lactone + H2O == D_Galactonate,
    D_galactonate + H == D_galactono_1_5_lactone + H2O,
    D_galactono_1_5_lactone + NADH + H == D_galactose + NAD,
    .6_phospho_D_galactono_1_5_lactone + NADH + H == D_galactose_6_phosphate + NAD,
    .6_phospho_D_gluconate + NAD == .6_phospho_2_dehydro_D_gluconate + NADH + H,
    .6_phospho_D_gluconate + NAD == D_ribulose_5_phosphate + CO2 + NADH,
    .6_oxohexanoate + NADH + H == .6_hydroxyhexanoate + NAD,
    CO2 + NADH == formate + NAD,
    D_alanine + NAD + H2O == pyruvate + NH4 + NADH + H,
    D_galactose + NAD == D_galactono_1_4_lactone + NADH + H,
    D_galactose + NADH + H == galactitol + NAD,
    D_glucarate + NADH + 2 * H == D_glucuronate + NAD + H2O,
    D_ribulose + NADH + H == D_arabinitol + NAD,
    D_ribulose + NADH + H == ribitol + NAD,
    D_ribulose_5_phosphate + CO2 + NADH == .6_phospho_D_gluconate + NAD,
    D_ribulose_5_phosphate + NADH + H == D_ribitol_5_phosphate + NAD,
    D_xylulose + NADH + H == D_arabinitol + NAD,
    D_xylulose + NADH + H == xylitol + NAD,
    D_xylose + NAD == D_xylono_1_5_lactone + NADH + H,
    FAD + NADH + 2 * H == FADH2 + NAD,
    FMNH2 + NAD == FMN + NADH + 2 * H,
    H2O2 + NADH + H == NAD + 2 * H2O,
    H2O2 + NAD == NADH + O2 + H,
    H2 + NAD == NADH + H,
    galactitol + NAD == D_galactose + NADH + H,
    galactitol + NAD == keto_D_tagatose + NADH + H,
    glycerol + NAD == dihydroxyacetone + NADH + H
);
let engine = vcell.model(vcell) |> engine.load(	
	inits = mass0(vcell, random = [10,50]),
	iterations       = 20, 
	time_resolutions = 3000, 	
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

bitmap(file = relative_work("large_mass.png"), width = 3600, height = 1920) {
    let p = ggplot(result_mass, padding = "padding: 5% 30% 10% 7%;");

    for(name in colnames(result_mass)) {
        if (name != "time.axis") {
            p = p + geom_line(aes(x = "time.axis", y = name), width = 8);
        }        
    }

    p;
}

bitmap(file = relative_work("large_flux.png"), width = 3600, height = 1920) {
    let p = ggplot(result_flux, padding = "padding: 5% 30% 10% 7%;");

    for(name in colnames(result_flux)) {
        if (name != "time.axis") {
            p = p + geom_line(aes(x = "time.axis", y = name), width = 8);
        }        
    }

    p;
}