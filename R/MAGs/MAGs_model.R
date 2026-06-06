imports "taxonomy_kit" from "metagenomics_kit";
imports "bifrost" from "seqtoolkit";

#' Compile the MAGs as virtual cell model
#' 
#' @param contigs a collection of the MAGs contigs assembly
#' @param tax_str taxonomy string of this MAGs contigs assembly data, example as:
#'    d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina sp001714685
#'
const MAGs_model = function(contigs, tax_str, outputdir = "./") {
    # make gene predictions from the MAGs contigs
    let genes = bifrost::prodigal(contigs, min.ORF.len = 90);

    if (is.character(tax_str)) {
        tax_str <- biom_string.parse(tax_str);
    }

    # export result to files
    write.csv(as.data.frame(result), file = file.path(outputdir, "gene_predicts.csv"));
    # export the gene prediction result to GFF3 format
    write.gff3(as.gff3(result), file = file.path(outputdir,"gene_predicts.gff3"));
    # export gene/protein fasta sequence to file
    write.fasta(as.genes(result), file = file.path(outputdir,"gene_predicts.fna"));
    write.fasta(as.proteins(result), file = file.path(outputdir,"protein_predicts.faa"));
}