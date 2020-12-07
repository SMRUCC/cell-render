imports ["geneExpression", "sampleInfo"] from "phenotype_kit";
imports "visualPlot" from "visualkit";
imports ["GSEA", "profiles"] from "gseakit";
imports ["dataset", "umap"] from "MLkit";

require(igraph);
require(igraph.render);

let umap_cluster_visual as function(matrix, outputdir) {

}

#' create clusters for biological function analysis
#' 
#' @param workspace the analysis workspace object, it should contains a 
#'     directory path list which can be used for read dep pvalue_cut 
#'     result.
#' 
let patterns_plot as function(workspace, matrix, outputdir) {
	matrix = getUnionDep(workspace, matrix); 



	lapply(workspace$analysis, compare_dir -> workspace :> create_pattern(compare_dir, outputdir));
}

#' create cluster for a specific analysis compare groups
#' 
#' @param workspace the workspace object that contains the directory path
#'     data for read pvalue_cut dep analysis result
#' @param compare_dir the analsysis design groups, is a list object that should contains
#'     a slot name \code{treatment} and a slot name \code{control}
#' 
let create_pattern as function(workspace, compare_dir, outputdir) {
	let pvalue_cut  = `${workspace$dirs$dep_analysis}/${as_label(compare_dir)}/pvalue_cut.csv`;
	let cluster_out = `${workspace$dirs$biological_analysis}/patterns/${as_label(compare_dir)}`;
	let kegg_background = read.background(`${outputdir}/annotation/kegg.Xml`);

	# removes all of the unnecessary information
	# create a expression matrix that contains
	# only of the protein expression normalization
	# data
	pvalue_cut = stripPvalue_cut(read.csv(pvalue_cut, row_names = 1));
				
	print(`[${as_label(compare_dir)}] previews of the different expression proteins:`);
	print(head(pvalue_cut));

	kegg_enrich(kegg_background, rownames(pvalue_cut), cluster_out);

	# run cmeans clustering
	let patterns = pvalue_cut 
	:> load.expr
	:> relative
	:> expression.cmeans_pattern(dim = [3,3], fuzzification = 2, threshold = 0.001)
	;
	
	# do data visualization
	patterns 
	:> plot.expression_patterns(size = [6000, 4500])
	:> save.graphics(`${cluster_out}/patterns.png`)
	;
	
	# and then save the cluster result matrix
	let matrix = patterns
	:> cmeans_matrix(kmeans_n = 5)
	;
	
	write.csv(matrix, file = `${cluster_out}/clusters.csv`);

	patterns = cluster.groups(matrix);

	str(patterns);

	for(groupKey in names(patterns)) {
		let members   = patterns[[groupKey]];
		let group_out = `${cluster_out}/cluster_${groupKey}`;

		print(head(members));

		# pvalue_cut[members, ]
		# :> dist 
		# :> btree(hclust = TRUE)
		# :> plot(       
		# 	size        = [3300, 30000], 
		# 	padding     = "padding: 200px 400px 200px 200px;", 
		# 	axis.format = "G2",
		# 	links       = "stroke: darkblue; stroke-width: 8px; stroke-dash: dash;",
		# 	pt.color    = "gray",
		# 	label       = "font-style: normal; font-size: 10; font-family: Bookman Old Style;",
		# 	ticks       = "font-style: normal; font-size: 12; font-family: Bookman Old Style;"
		# )
		# :> save.graphics(`${group_out}/deps.png`)
		# ;

		kegg_enrich(kegg_background, members, group_out);
		writeLines(members, con = `${group_out}/deps.txt`);
	}
}

let kegg_enrich as function(kegg_background, members, group_out) {
	# run GSEA
	let enrich = kegg_background 
	:> enrichment(members, showProgress = FALSE) 
	:> enrichment.FDR;

	enrich :> write.enrichment(
		file   = `${group_out}/kegg.csv`, 
		format = "KOBAS"
	);

	enrich 
	:> as.KOBAS_terms 
	:> KEGG.enrichment.profile(top = 10)
	:> category_profiles.plot(
		size       = [2700,2400], 
		dpi        = 100, 
		title      = "KEGG enrichment", 
		axis_title = "-log10(pvalue)"
	)
	:> save.graphics(file = `${group_out}/kegg.png`)
	;
}