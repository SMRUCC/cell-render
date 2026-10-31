imports "bioseq.blast" from "seqtoolkit";

#' run diamond blastp search via the GCModeller managed method
const clr_diamond = function(q, s) {
    diamond("Fast") |> diamond_search(q,s);
}