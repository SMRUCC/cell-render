imports "visualPlot" from "visualkit";
imports ["geneExpression", "sampleInfo"] from "phenotype_kit";

let run_dep as function(matrix, sampleInfo, compares, output_dir) {
	let as_label = compare_dir -> `${compare_dir$treatment} vs ${compare_dir$control}`;
	let	sampleinfo = read.sampleinfo(sampleInfo); 

	for(compare_dir in compares) {
		print(`run dep analysis of '${as_label(compare_dir)}'...`);
	
		let dep = deg.t.test(matrix, sampleinfo, compare_dir$treatment, compare_dir$control);
		let compare_out = `${output_dir}/analysis/03.dep_analysis/${as_label(compare_dir)}`;
		
		write.csv(dep, file = `${compare_out}/pvalue.csv`);
		volcano.plot(dep, 
			size = "1400,1600", 
			title = `Volcano plot of ${as_label(compare_dir)}`) 
		:> save.graphics(file = `${compare_out}/volcano.png`)
		;
		
		print("done!");
	}
}