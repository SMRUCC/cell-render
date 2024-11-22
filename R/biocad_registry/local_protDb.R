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
const local_protDb = function(cad_registry, dbfile) {
    let prot_term = 3;
    let page_size = 25000;
    let offset = 1;
    let stream = open.fasta(dbfile, read = FALSE);
    let total_size = cad_registry |> table("molecule") 
    |> where(type = prot_term) 
    |> count();
    let total_pages = as.integer(total_size / page_size) + 1;

    for(let page in 1:total_pages) {
        print(`load data page ${page}:`);

        offset <- (page - 1) * page_size;
        page <- cad_registry 
        |> table("molecule") 
        |> left_join("sequence_graph")
        |> on(sequence_graph.molecule_id = molecule.id)
        |> left_join("db_xrefs")
        |> on(db_xrefs.obj_id = molecule.id, db_xrefs.db_key in [77 , 79])
        |> where(molecule.type = prot_term)
        |> group_by("molecule.id")
        |> limit(offset, page_size)
        |> select(
            "`molecule`.id",
            "MIN(sequence) AS prot_seq",
            "MIN(note) AS prot",
            "GROUP_CONCAT(DISTINCT xref) AS ec_number")
        ;

        # print(page);
        # stop();

        page[, "ec_number"] <- ifelse(nchar(page$ec_number) == 0, "-", page$ec_number);
        page[, "ref"] <- `${page$id} ${page$ec_number}`;
        page <- as.list(page, byrow = TRUE);
        page 
        |> tqdm()
        |> sapply(i -> fasta(i$prot_seq, attrs = [i$ref, i$prot])) 
        |> write.fasta(file = stream, filter_empty = TRUE)
        ;
    }

    close(stream);
}