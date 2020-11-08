imports "../package/01_background.R";

setwd(dirname(!script$dir));

let output_dir as string = "demo";

makePtf(`${output_dir}/annotation/uniprot-taxonomy_3702.xml`, `${output_dir}/annotation/background.ptf`);