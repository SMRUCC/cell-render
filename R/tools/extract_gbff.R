imports "GenBank" from "seqtoolkit";
imports "bioseq.fasta" from "seqtoolkit";
imports "annotation.genomics" from "seqtoolkit";

#' Extract Upstream Sequences and Annotations from a GenBank File
#'
#' This is a helper function that parses a GenBank assembly file to extract 
#' genomic sequences, gene features, and the upstream regions relative to the 
#' Transcription Start Sites (TSS). It is designed to prepare data for 
#' downstream transcript regulation network analysis.
#'
#' The function performs extensive file I/O operations. It extracts the raw 
#' genomic FASTA, calculates the TSS upstream loci for each gene, cuts these 
#' loci from the genome, and exports both the sequences and tabular context 
#' files to the specified working directory.
#'
#' @param src A character string specifying the file path to the input 
#'   GenBank (`.gb`, `.gbff`, or `.gbk`) file.
#' @param workdir A character string specifying the directory path where 
#'   all output files will be saved. Defaults to `"./"` (current working directory).
#' @param upstream_size An integer specifying the length of the upstream region 
#'   (in base pairs) to extract relative to the gene TSS. Defaults to `150`.
#' @param tag_genbank_accid A logical value. If `TRUE`, the GenBank accession 
#'   ID is appended as the source tag in the FASTA headers of the extracted 
#'   upstream sequences. If `FALSE` (default), the original sequence headers 
#'   are used as the source tag.
#' @param verbose A logical value. If `TRUE` (default), progress messages 
#'   and the titles of the extracted upstream loci are printed to the console.
#'
#' @return Returns `invisible(NULL)`. This function is called primarily for 
#'   its side effects: generating and saving output files to the `workdir`.
#'
#' @section Output Files:
#' The function generates the following files in the `workdir`:
#' \itemize{
#'   \item{\code{genes.csv}: }{A CSV file containing the extracted gene feature table.}
#'   \item{\code{context.txt}: }{A PTT (Protein Table) format file containing the genomic context.}
#'   \item{\code{source.fasta}: }{A FASTA file of the complete raw genomic sequence.}
#'   \item{\code{upstream_locis.fasta}: }{A FASTA file containing the extracted upstream 
#'     sequences. Headers are formatted as `[Gene_ID] [Source_Tag]`.}
#' }
#'
#' @importFrom seqtoolkit GenBank
#' @importFrom seqtoolkit bioseq.fasta
#' @importFrom seqtoolkit annotation.genomics
#' @export
const extract_gbff = function(src, workdir = "./", 
                              upstream_size = 150, 
                              tag_genbank_accid = FALSE, 
                              verbose = TRUE) {

    # load genbank assembly file from a given file path
    let gbk = read.genbank(src);
    # extract the raw genomics fasta sequence
    let genomics_seq = GenBank::origin_fasta(gbk);
    # extract the gene features from the genbank assembly object
    let genes = genome.genes(genome = gbk);
    let gene_ids = [genes]::Synonym;
    let genbank_id = gbk |> GenBank::accession_id();

    print(`target genome genbank accession id: ${genbank_id}.`);
    print("get genes table:");
    print(gene_ids);
    print("bp size for parse the gene upstream loci:");
    str(upstream_size);

    let locis = genes 
    # extract gene TSS upstream region and then assign the list name with gene ids
    |> upstream(length = upstream_size, is_relative_offset = TRUE) 
    |> lapply(l => l, names = gene_ids)
    |> tqdm()
    # cast each gene TSS upstream location as site fasta sequence
    # by cut site sequence from genomics sequence via the
    # given TSS upstream location data
    |> lapply(function(loci, i) {
        let fa = cut_seq.linear(genomics_seq, loci, nt_auto_reverse = TRUE);
        let id = gene_ids[i];
        let source_tag = {
            if (tag_genbank_accid) {
                genbank_id;
            } else {
                [fa]::Headers;
            }
        }

        # tag the corresponding gene id to the
        # loci site headers
        fasta.headers(fa) <- append(id, source_tag);
        fa;
    })
    ;

    if (verbose) {
        print("view upstream locis:");
        print(fasta.titles(locis));
    }

    # export genomics context elements as feature table.
    write.csv(genes, file = `${workdir}/genes.csv`);
    # export PTT genomics context tabular file
    write.PTT_tabular(gbk, file = `${workdir}/context.txt`);

    # extract sequence data and the gene context data for the
    # downstream transcript regulation network analysis
    write.fasta(genomics_seq, file = `${workdir}/source.fasta`);
    write.fasta(locis, file = `${workdir}/upstream_locis.fasta`);

    invisible(NULL);
}