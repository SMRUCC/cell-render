imports "visualPlot" from "visualkit";
imports "geneExpression" from "phenotype_kit";

let run_dep as function(matrix, sampleinfo, compares, output_dir) {
	let as_label = compare_dir -> `${compare_dir$treatment} vs ${compare_dir$control}`;

	for(compare_dir in compares) {
		let dep = deg.t.test(matrix, sampleinfo, compare_dir$treatment, compare_dir$control);
		let compare_out = `${output_dir}/${as_label(compare_dir)}`;
		
		write.csv(dep, file = `${compare_out}/pvalue.csv`);
		volcano.plot(dep) :> save.graphics(file = `${compare_out}/volcano.png`);
	}
}