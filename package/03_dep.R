imports "visualPlot" from "visualkit";
imports ["geneExpression", "sampleInfo"] from "phenotype_kit";

let as_label = compare_dir -> `${compare_dir$treatment} vs ${compare_dir$control}`;

let run_dep as function(workspace, matrix) {
	let	sampleinfo = workspace$sample_info; 

	for(compare_dir in workspace$analysis) {
		print(`run dep analysis of '${as_label(compare_dir)}'...`);
	
		let dep = deg.t.test(matrix, sampleinfo, compare_dir$treatment, compare_dir$control, level = 1.25, FDR = 1);
		let compare_out = `${workspace$dirs$dep_analysis}/${as_label(compare_dir)}`;
		let pvalue_cut = as.data.frame(dep[sapply(dep, prot -> as.object(prot)$isDEP)]);
		
		rownames(pvalue_cut) = make.names(rownames(pvalue_cut), unique = TRUE, allow_ = TRUE);
		
		write.csv(dep, file = `${compare_out}/pvalue.csv`);
		write.csv(pvalue_cut, file = `${compare_out}/pvalue_cut.csv`);
		
		# data visualization
		volcano.plot(dep, 
			size = "1400,1600", 
			title = `Volcano plot of ${as_label(compare_dir)}`) 
		:> save.graphics(file = `${compare_out}/volcano.png`)
		;
		
		print("done!");
	}
}