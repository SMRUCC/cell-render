require(GCModeller);
require(umap);
require(igraph);

imports ["annotation.workflow", "annotation.terms"] from "seqtoolkit";
imports "taxonomy_kit" from "metagenomics_kit";
imports "OTU_table" from "metagenomics_kit";

const diamond_embedding = function(diamond_result, workdir = "./", union_contigs = 250) {
    let rawdata = read_m8(diamond_result, stream = TRUE);
    let ec_terms = file.path(workdir, "ec_terms.csv");
    let metabolic_data  = file.path(workdir, "genomes_metabolic.jsonl");
    let ec_emebedding = file.path(workdir, "metabolic_embedding.csv");
    let stream = open.stream(ec_terms, type = "terms");

    stream.flush(m8_metabolic_terms(rawdata), stream);
    stream = open.stream(ec_terms,type = "terms", ioRead = TRUE);

    let genomes = make_vectors(stream, stream = TRUE);

    write_genomes_jsonl(genomes, file = metabolic_data);

    let models = read.jsonl(file = metabolic_data, what = "genome_vector");
    let vec = models |> tfidf_vectorizer(union_contigs = union_contigs, hierarchical  = TRUE);

    write.csv(vec, file = ec_emebedding );
}


