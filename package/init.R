const summary_info as string = "
raw data summary

this folder contains the raw data summary analysis, 
includes raw data overviews plot, GO and KEGG 
annotations.
";

const stat_info as string = "
statistical analysis

this folder contains the statistical analysis of your
data samples, mainly includes the PCA analysis, etc.
";

const dep_info as string = "
different expression proteins

this folder contains the different expression proteins 
analysis result between your sample groups.
";

const biological_info as string = "
biological function analysis

this folder contains the biological function analysis of
your different expression proteins, mainly contains of the
GO and KEGG enrichment analysis result.
";

const ppi_info as string = "
Protein-Protein interaction analysis

this folder contains the result of the PPI analysis, mainly
includes of the different expression protein correlation 
calculation analysis result of PCC and SPCC. And also includes
a correlation network visualization in this analysis.
";

const folders = list(
	"01.summary"             = summary_info, 
	"02.stat"                = stat_info, 
	"03.dep_analysis"        = dep_info, 
	"04.biological_analysis" = biological_info, 
	"05.ppi_analysis"        = ppi_info
);

let init_workspace as function(output_dir as string) {
	let workspace = list();

	for(name in names(folders)) {
		let dir as string = `${output_dir}/analysis/${name}`;
		let workspace_name = strsplit(name, ".", fixed = TRUE)[2];
		
		workspace[[workspace_name]] = dir;
		
		dir.create(dir);
		
		# write readme text file
		folders[[name]] 
		:> trim("\n") 
		:> writeLines(con = `${dir}/readme.txt`)
		;
	}
	
	list(dirs = workspace);
}