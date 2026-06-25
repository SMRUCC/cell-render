imports "metaTraits" from "metagenomics_kit";

#' run workflow of metaTraits for make phenotype traits annotation
#' 
const metaTraits = function(diamond) {
    let pfam = read_m8(diamond);
    let traitar = metaTraits::load.trait_models(repo = `${@datadir}/phypat+PGL/`);
    let votes = metaTraits::make_predicts(pfam, bit_score = 25, evalue = 0.01);
    let pheno_traits = metaTraits::phenotype_result(votes, models = traitar );

    return(pheno_traits);
}