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
#' @param src A character string specifying the file path to the GenBank
#'   source file.
#' @param workdir A character string specifying the working directory where
#'   all extracted files will be saved. Defaults to \code{"./"}.
#' @param upstream_size An integer specifying the number of base pairs
#'   upstream of the TSS to extract. Defaults to 300.
#' @param verbose Logical. If \code{TRUE}, prints additional debug
#'   information including the upstream locus titles. Defaults to \code{FALSE}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of writing the following files to \code{workdir}:
#'   \describe{
#'     \item{\code{genes.csv}}{Gene feature annotation table extracted from
#'       the GenBank file.}
#'     \item{\code{context.txt}}{PTT-format genomic context tabular file.}
#'     \item{\code{source.fasta}}{Raw genomic nucleotide FASTA sequences.}
#'     \item{\code{upstream_locis.fasta}}{TSS upstream region sequences
#'       with gene locus tag headers.}
#'   }
#'
#' @details
#' The extraction pipeline:
#' \enumerate{
#'   \item Reads the GenBank file and extracts the full genomic sequence.
#'   \item Parses gene features and their coordinates.
#'   \item Calculates the TSS position for each gene based on strand orientation.
#'   \item Extracts the upstream region of the specified size from each TSS.
#'   \item Tags each upstream sequence with its corresponding gene locus tag.
#'   \item Exports all data to the working directory.
#' }
#'
#' @seealso \code{\link{make_genbank_proj}} for the higher-level workflow
#'   that calls this function, \code{\link{tfbs_motif_scanning}} which
#'   consumes the \code{upstream_locis.fasta} output.
#'
#' @examples
#' \dontrun{
#' extract_gbff(
#'   src = "data/sequence.gbff",
#'   workdir = "results/extracted",
#'   upstream_size = 500,
#'   verbose = TRUE
#' )
#' }
#'
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