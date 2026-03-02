require(GCModeller);

imports ["annotation.workflow", "annotation.terms"] from "seqtoolkit";

let rawdata = read_m8(relative_work("ec_numbers_results.tsv"), stream = TRUE);
let stream = open.stream(relative_work("ec_terms.csv"),type = "terms");

stream.flush(m8_metabolic_terms(rawdata), stream);