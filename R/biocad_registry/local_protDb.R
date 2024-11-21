#' Create local fasta database for protein alignment
#' 
#' @details this function will pull all protein sequence from the biocad_registry, 
#' and then build local fasta sequence database with header formats:
#' 
#' ```
#' > cad_id ec_number|function description 
#' ```
#' 
#' a placeholder symbol ``-`` will be placed after the cad_id if the ``ec_number`` 
#' of the protein is missing. for multiple ec_number, a comma seperator will be 
#' used for join the ec_number list.
#' 
const local_protDb = function(cad_registry) {
    let page_size = 1000;
    let offset = 1;
    let 

    for(let page in tqdm(1:10000)) {
        offset <- (page - 1) * page_size;
        page <- cad_registry 
        |> table("molecule") 
        |> left_join("sequence_graph")
        |> on(sequence_graph.molecule_id = molecule.id)
        |> left_join("db_xrefs")
        |> on(db_xrefs.obj_id = molecule.id, db_key in [77 , 79])
        |> where(molecule.type = 3)
        |> group_by("molecule.id")
        |> limit(offset, page_size)
        |> select(
            "`molecule`.id",
            "MIN(sequence) AS prot_seq",
            "MIN(note) AS prot",
            "GROUP_CONCAT(DISTINCT xref) AS ec_number")
        ;

        print(page);
        stop();
    }
}