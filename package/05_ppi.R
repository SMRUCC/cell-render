let dep_correlations as function(workspace, matrix, outputdir) {  
    let union_dep = lapply(workspace$analysis, compare_dir -> getDepGeneId(workspace, compare_dir));

    union_dep = unique(unlist(union_dep));
    output_dir = `${output_dir}/${workspace$dirs$ppi_analysis}`;

    print(`there are ${length(union_dep)} union dep from all of your analysis groups:`);
    print(head(union_dep));

    matrix = as.data.frame(matrix)[union_dep, ];

    print("previews of your union dep raw expression data matrix:");
    print(head(matrix));

}

let getDepGeneId as function(workspace, compare_dir) {
    let pvalue_cut as string = `${workspace$dirs$dep_analysis}/${as_label(compare_dir)}/pvalue_cut.csv`;
    let data = read.csv(pvalue_cut, row_names = 1);

    rownames(data);
}