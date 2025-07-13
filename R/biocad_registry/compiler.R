#' Compile the virtual cell model from ncbi genbank file
#' 
#' @param cad_registry the mysql connection to the biocad registry
#' @param gbff the file path to the template ncbi genbank file
#' 
const compile_genbank = function(cad_registry, gbff) {
    imports "GenBank" from "seqtoolkit";

    if (is.character(gbff)) {
        gbff <- read.genbank(gbff);
    }

    let ncbi_taxid = taxon_id(gbff);
    let cellular_id = accession_id(gbff);
    let taxon = taxonomy_lineage(gbff);
    
    cad_registry |> Builder::create_modelfile(
        genes = gbff |> GenBank::as_tabular(ORF = FALSE),
        taxid = ncbi_taxid,
        cellular_id = cellular_id,
        taxonomy = taxon
    );
}