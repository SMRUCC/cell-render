require(GCModeller);
require(umap);

imports "taxonomy_kit" from "metagenomics_kit";
imports "annotation.terms" from "seqtoolkit";

let models = read.jsonl(file = relative_work("genomes_metabolic.jsonl"), what = "genome_vector");
let vec = models |> tfidf_vectorizer(union_contigs = 250, hierarchical  = TRUE);

write.csv(vec, file = relative_work("metabolic_embedding.csv"));