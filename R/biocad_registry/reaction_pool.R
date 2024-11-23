#' pull reaction data from cad registry to local
#' 
const reaction_pool = function(cad_registry) {
    let total = cad_registry |> table("reaction") |> count();
    let page_size = 1000;
    let total_pages = as.integer(total / page_size) + 1;
    let offset = 0;

    for(let page in 1:total_pages) {
        print(`processing of data page: ${page}/${page_size}...`);

        offset <- (page - 1) * page_size;
        page <- cad_registry |> table("reaction") 
            |> limit(offset, page_size) 
            |> select()
            ;

        for(let r in as.list(page, byrow = TRUE) |> tqdm()) {
            page <- cad_registry |> table("reaction_graph") 
                |> left_join("sequence_graph")
                |> on(sequence_graph.molecule_id = reaction_graph.molecule_id)
                |> where(reaction = r$id)
                |> select(
                    sequence_graph.molecule_id,
                    db_xref,
                    role,
                    "sequence AS smiles",
                    embedding
                )
                ;

            if (any(nchar(page$embedding) == 0)) {
                next;
            }

            print(page);
            stop();
        }
    }
}