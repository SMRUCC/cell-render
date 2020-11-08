imports "../package/01_background.R";

setwd(dirname(!script$dir));

let output_dir as string = "demo";
let overrides as boolean = ?"--overrides" || FALSE;

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