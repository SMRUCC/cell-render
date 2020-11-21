imports "visualPlot" from "visualkit";
imports ["geneExpression", "sampleInfo"] from "phenotype_kit";

#' convert the analysis compares groups as the 
#' folder name label
#' 
#' @param compare_dir the analysis groups data, is a list object that should
#'    contains a slot named \code{treatment} and a slot named \code{control}.
#'    both slot value should be the \code{sample_info} value in the analysis
#'    samples.
#' 
let as_label = compare_dir -> `${compare_dir$treatment} vs ${compare_dir$control}`;

#' run dep analysis for the raw data base on each
#' given analysis compare groups 
#' 
#' @param matrix the protein expression data, it should be normalized by
#'    total sum of the proteins' peak area
#' @param workspace a analysis workspace object that contains the directory 
#'    for save data and also contains the analysis parameters, like threshold
#'    values, sample information and compare analaysis designs. 
#' 
let run_dep as function(workspace, matrix) {
	let	sampleinfo = workspace$sample_info; 

	for(compare_dir in workspace$analysis) {
		workspace :> dep_calc(matrix, compare_dir, sampleinfo);
	}
}

#' dep calculation
#' 
#' @param workspace an analysis workspace object that contains directory path
#'   for save result data, andalso it should contains the analysis parameters
#' @param matrix the normalized protein expression matrix object
#' @param compare_dir the analysis compres group
#' @param sampleinfo the sample information of the sample columns in the protein 
#'    expression \code{matrix} data.
#' 
let dep_calc as function(workspace, matrix, compare_dir, sampleinfo) {
	let log2FC_level as double = workspace$args$log2FC_level || 1.25;
	let FDR_threshold as double = workspace$args$FDR || 1;

	print(`run dep analysis of '${as_label(compare_dir)}'...`);

	# apply of the t.test for the dep analysis
	let dep = deg.t.test(
		matrix, sampleinfo, compare_dir$treatment, compare_dir$control, 
		level = log2FC_level, 
		FDR = FDR_threshold
	);
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