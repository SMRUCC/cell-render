imports ["geneExpression", "sampleInfo"] from "phenotype_kit";
imports "visualPlot" from "visualkit";

#' create clusters for biological function analysis
#' 
#' @param workspace the analysis workspace object, it should contains a 
#'     directory path list which can be used for read dep pvalue_cut 
#'     result.
#' 
let patterns_plot as function(workspace) {
	lapply(workspace$analysis, compare_dir -> workspace :> create_pattern(compare_dir));
}

#' create cluster for a specific analysis compare groups
#' 
#' @param workspace the workspace object that contains the directory path
#'     data for read pvalue_cut dep analysis result
#' @param compare_dir the analsysis design groups, is a list object that should contains
#'     a slot name \code{treatment} and a slot name \code{control}
#' 
let create_pattern as function(workspace, compare_dir) {
	let pvalue_cut  = `${workspace$dirs$dep_analysis}/${as_label(compare_dir)}/pvalue_cut.csv`;
	let cluster_out = `${workspace$dirs$biological_analysis}/${as_label(compare_dir)}`;
	
	# removes all of the unnecessary information
	# create a expression matrix that contains
	# only of the protein expression normalization
	# data
	pvalue_cut = stripPvalue_cut(read.csv(pvalue_cut, row_names = 1));
				
	print(`[${as_label(compare_dir)}] previews of the different expression proteins:`);
	print(head(pvalue_cut));
	
	pvalue_cut 
	:> dist 
    :> btree(hclust = TRUE)
    :> plot(       
        size        = [3300, 23000], 
        padding     = "padding: 200px 400px 200px 200px;", 
        axis.format = "G2",
        links       = "stroke: darkblue; stroke-width: 8px; stroke-dash: dash;",
        pt.color    = "gray",
        label       = "font-style: normal; font-size: 10; font-family: Bookman Old Style;",
        ticks       = "font-style: normal; font-size: 12; font-family: Bookman Old Style;"
    )
    :> save.graphics(`${cluster_out}/deps.png`)
    ;

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
	patterns
	:> cmeans_matrix(kmeans_n = 3)
	:> write.csv(file = `${cluster_out}/clusters.csv`)
	;
}