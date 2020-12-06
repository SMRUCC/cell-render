# run init of the workspace
imports "../package/init.R";
# data analysis and pre-processing components
imports "../package/01_background.R";
imports "../package/01_raw.R";
imports "../package/03_dep.R";
imports "../package/04_biological.R";
imports "../package/05_ppi.R";

setwd(dirname(!script$dir));

let output_dir as string = "demo";
let overrides as boolean = ?"--overrides" || FALSE;
let sample_info as string = `${output_dir}/sampleInfo.csv`;
let goDb as string = ?"--go" || "P:/go.obo";

sink(file = `${output_dir}/analysis/pipeline.log`);

# create workspace folders
let workspace = init_workspace(output_dir);

workspace$args = list(log2FC_level = 2);
workspace$analysis = list(
	a = list(treatment = "C6", control = "C9"),
	b = list(treatment = "C6", control = "I56"),
	c = list(treatment = "C6", control = "I59"),
	d = list(treatment = "C6", control = "I86"),
	e = list(treatment = "C6", control = "I89")
);
workspace$sample_info = read.sampleinfo(sample_info); 

print("view of the analysis workspace object in current pipeline:");
str(workspace);

let background_ptf as string = `${output_dir}/annotation/background.ptf`;
let HTS as string = `${output_dir}/raw/uniprot.csv`;
let uniprot_src = `${output_dir}/annotation/uniprot-taxonomy_3702.xml`;

if (overrides) {
	print(`annotation background '${background_ptf}' will be overrides!`);
}

# save background annotation data
if (overrides || !file.exists(background_ptf)) {
	# makePtf(uniprot_src, background_ptf);

	# unify the gene id to uniprot protein id.
	`${output_dir}/raw/all_counts.csv`
	:> unifyId(background_ptf) 
	:> write.csv(file = HTS)
	;

	GSEAbackground(background_ptf, `${output_dir}/annotation`); 
}

let annotations = read.csv(HTS, row_names = 1)
:> protein_annotations(ptf = background_ptf)
;

# stage 01, raw sample data analysis
# includes protein function annotations in current expression samples
annotations
:> write.csv(file = `${workspace$dirs$summary}/protein.annotations.csv`)
;
annotations
:> go_summary(goDb, `${workspace$dirs$summary}/GO`)
;

workspace :> hist_samples(matrix = read.csv(HTS, row_names = 1));
# run dep analysis and data visualization of the dep
workspace :> run_dep(matrix = load.expr(HTS, rm_ZERO = TRUE));
# create cluster for biological function analysis
workspace :> patterns_plot(output_dir);
workspace :> dep_correlations(matrix = load.expr(HTS, rm_ZERO = TRUE), output_dir = output_dir);

print("Workflow finished!");

sink();