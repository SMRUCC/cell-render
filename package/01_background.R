imports "http" from "webKit";
imports ["ptfKit", "labelfree", "summary"] from "proteomics_toolkit";
imports "uniprot" from "seqtoolkit";
imports "geneExpression" from "phenotype_kit";
imports "file" from "gokit";
imports "background" from "gseakit";
imports "repository" from "kegg_kit";

# script for create background annotation data files
# the annotation data is request from the uniprot database
# via web api
const url_template as string = "https://www.uniprot.org/uniprot/?query=taxonomy:%s&format=xml&force=true&compress=yes";
const keggMaps as string = "E:\biodeep\biodeepdb_v3\KEGG\br08901_pathwayMaps";

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
	:> uniprot.ptf(
		keys = ["KEGG","KO","GO","Pfam","RefSeq","EC","InterPro","BioCyc","eggNOG","EMBL","STRING","EnsemblPlants","Proteomes", "Araport"], 
		includesNCBITaxonomy = TRUE,
		scientificName = TRUE
	)
	:> save.ptf(file = save)
	;
}

#' unify the various annotation id into the uniprot id
#'
#' @param raw the matrix of the labelfree result data.
#'
#' @returns a raw matrix with all protein id unify as the uniprot id.
#'
let unifyId as function(raw, ptf) {
	let genes = raw 
	:> load.expr(rm_ZERO = TRUE) 
	:> as.generic 
	# apply of the total sum normalization 
	# for the label free samples data
	:> sample.normalize
	;
	
	ptf 
	:> load.ptf 
	:> as.uniprot_id(genes)
	;
}

#' exports of the protein annotation table data
#' 
#' @param raw the expression sample data, the protein id will be used for 
#'     pick of the protein annotation data in the unify annotation model
#' 
#' @param ptf the GCModeller unify annotation data model
#' 
let protein_annotations as function(raw, ptf) {
	let geneIDs as string = rownames(raw);

	print("previews part of the unify protein ids:");
	print(head(geneIDs));

	ptf 
	:> load.ptf 
	:> protein.annotations(geneIDs)
	;
}

#' GO annotation summary of the proteins in current sample data
#'
#' @param annotations a protein annotation table data that created by 
#'                    \code{protein_annotations}
#' @param goDb the file path of the GO obo database file
#' @param outputdir the directory for save the count table and the bar 
#'                  plot image file.
#'
let go_summary as function(annotations, goDb, outputdir) {
	let profiles = annotations 
	:> proteins.GO(goDb = read.go_obo(goDb))
	;

	# namespace go_term description counts
	as.data.frame(profiles, type = "go") 
	:> write.csv(file = `${outputdir}/counts.csv`, row_names = FALSE)
	;

	profiles
	:> profileSelector(selects = "desc:13")
	:> category_profiles.plot(
		title = "GO profiles", 
		axis_title = "Number Of Proteins",
		dpi = 150,
		size = [2300, 2200]
	)
	:> save.graphics(file = `${outputdir}/plot.png`)
	;
}

let GSEAbackground as function(ptf, outputdir) {
	let annotations = ptf :> load.ptf; 
	let kegg_maps = load.maps(keggMaps, rawMaps = FALSE);
	let kegg_xml as string = `${outputdir}/kegg.Xml`;

	print("kegg background model will be saved at location:");
	print(kegg_xml);

	annotations 
	:> KO.background(kegg_maps, size = length(as.object(annotations)$proteins))
	:> write.background(file = kegg_xml)
	;
}