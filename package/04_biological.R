imports ["geneExpression", "sampleInfo"] from "phenotype_kit";
imports "visualPlot" from "visualkit";

let patterns_plot as function(workspace) {
	for(compare_dir in compares) {
		let pvalue_cut  = `${workspace$dirs$dep_analysis}/${as_label(compare_dir)}/pvalue_cut.csv`;
		let cluster_out = `${workspace$dirs$biological_analysis}/${as_label(compare_dir)}`;
		
		pvalue_cut = read.csv(pvalue_cut);
		pvalue_cut[, "FC.avg"]  = NULL;
		pvalue_cut[, "p.value"] = NULL;
		pvalue_cut[, "is.DEP"]  = NULL;
		pvalue_cut[, "log2FC"]  = NULL;
		pvalue_cut[, "FDR"]     = NULL;
		
		let patterns = pvalue_cut 
		:> load.expr
		:> relative
		:> expression.cmeans_pattern(dim = [3,3], fuzzification = 2, threshold = 0.001)
		;
		
		patterns 
		:> plot.expression_patterns()
		:> save.graphics(`${cluster_out}/patterns.png`)
		;
		
		patterns
		:> cmeans_matrix(kmeans_n = 3)
		:> write.csv(file = `${cluster_out}/clusters.csv`)
		;
	}	
}
