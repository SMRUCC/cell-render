require(GCModeller);

imports ["annotation.workflow", "annotation.terms"] from "seqtoolkit";

let stream = open.stream(relative_work("ec_terms.csv"),type = "terms", ioRead = TRUE);
let genomes = make_vectors(stream, stream = TRUE);

write_genomes_jsonl(genomes, file = relative_work("genomes_metabolic.jsonl"));
