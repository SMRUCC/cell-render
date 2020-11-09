const summary_info as string = "
raw data summary

this folder contains the raw data summary analysis, 
includes raw data overviews plot, GO and KEGG 
annotations.
";

const stat_info as string = "

";

const dep_info as string = "

";

const biological_info as string = "

";

const ppi_info as string = "

";

const folders = list(
	"01.summary"             = summary_info, 
	"02.stat"                = stat_info, 
	"03.dep_analysis"        = dep_info, 
	"04.biological_analysis" = biological_info, 
	"05.ppi_analysis"        = ppi_info
);

let init_workspace as function(output_dir as string) {
	for(name in names(folders)) {
		let dir as string = `${output_dir}/analysis/${name}`;
		
		dir.create(dir);
		
		# write readme text file
		folders[[name]] 
		:> trim("\n") 
		:> writeLines(con = `${dir}/readme.txt`)
		;
	}
}