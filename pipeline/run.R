# run init of the workspace
imports "../package/init.R";
# data analysis and pre-processing components
imports "../package/01_background.R";
imports "../package/03_dep.R";

setwd(dirname(!script$dir));

let output_dir as string = "demo";
let overrides as boolean = ?"--overrides" || FALSE;
let sample_info as string = `${output_dir}/sampleInfo.csv`;
let analysis = list(
	a = list(treatment = "C6", control = "C9"),
	b = list(treatment = "C6", control = "I56"),
	c = list(treatment = "C6", control = "I59"),
	d = list(treatment = "C6", control = "I86"),
	e = list(treatment = "C6", control = "I89")
);

let background_ptf as string = `${output_dir}/annotation/background.ptf`;
let HTS as string = `${output_dir}/raw/uniprot.csv`;

if (overrides) {
	print(`annotation background '${background_ptf}' will be overrides!`);
}

if (overrides || !file.exists(background_ptf)) {
	makePtf(`${output_dir}/annotation/uniprot-taxonomy_3702.xml`, background_ptf);
}

`${output_dir}/raw/all_counts.csv`
:> unifyId(background_ptf) 
:> write.csv(file = HTS)
;

init_workspace(output_dir);

run_dep(
	matrix     = load.expr(HTS), 
	sampleInfo = sample_info, 
	compares   = analysis, 
	output_dir = output_dir
);