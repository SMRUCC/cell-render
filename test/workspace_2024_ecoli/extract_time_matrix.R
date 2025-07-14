require(GCModeller);

imports "rawXML" from "vcellkit";
imports "geneExpression" from "phenotype_kit";

setwd(@dir);

let rawdata = open.vcellPack(file = "./result.vcellPack", mode = "read");

for(let term in ["Transfer-RNA"  "Ribosomal-RNA" "Message-RNA" "Component-RNA" "Polypeptide" "Protein" "Metabolite"]) {
    let metabolome = time.frames(rawdata, module=term, symbol_name= TRUE);

    metabolome |> write.expr_matrix(file =`./vcell_sim/${term}.csv`, id = "molecules");

    require(ggplot);

    metabolome = as.data.frame(metabolome);
    metabolome[,"time.axis"] = 1:nrow(metabolome); 

    bitmap(file = `./vcell_sim/${term}.png`, width = 3600, height = 1920) {
        let p = ggplot(metabolome, padding = "padding: 5% 30% 10% 7%;");

        for(name in colnames(metabolome)) {
            if (name != "time.axis") {
                p = p + geom_line(aes(x = "time.axis", y = name), width = 8);
            }        
        }

        p;
    }
}

