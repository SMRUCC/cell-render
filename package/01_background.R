imports "http" from "webKit";
imports "ptfKit" from "proteomics_toolkit";
imports "uniprot" from "seqtoolkit";

# script for create background annotation data files
# the annotation data is request from the uniprot database
# via web api
const url_template as string = "https://www.uniprot.org/uniprot/?query=taxonomy:%s&format=xml&force=true&compress=yes";

let request_uniprot as function(taxid as string, save as string) {
	let url as string = sprintf(url_template, taxid);
	
	print("will downloads annotation background dataset from uniprot:");
	print(url);
	
	wget(url, save);
}

#' Make annotation from uniprot xml database files
#'
#' @param uniprot the file path of the uniprot xml database file
#'
let makePtf as function(uniprot as string, save as string) {
	uniprot
	:> open.uniprot
	:> uniprot.ptf(keys = ["KEGG","KO","GO","Pfam","RefSeq","EC","InterPro","BioCyc","eggNOG","EMBL","STRING","EnsemblPlants","Proteomes"])
	:> save.ptf(file = save)
	;
}

#' unify the various annotation id into the uniprot id
#'
#' @param raw the matrix of the labelfree result data.
#'
#' @returns a raw matrix with all protein id unify as the uniprot id.
#'
let unifyId as function(raw) {

}