require(CellRender);

let repo = "K:\bacterias\ncbi_dataset\data";

for(let file in list.files(repo, pattern = "*.gbff", recursive = TRUE)) {
    CellRender::extract_gbff(file, 
        workdir = dirname(file), 
        upstream_size = 150, 
        verbose = TRUE);
}

