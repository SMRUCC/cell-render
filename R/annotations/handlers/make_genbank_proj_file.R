#' Initialize a Single GenBank Project and Export Sequence Files
#'
#' Core processing function that ingests a GenBank file, wraps it into a 
#' standardized \code{project} object via \code{CellRender}, saves the 
#' project state to disk, and extracts essential FASTA files required for 
#' downstream bioinformatics analyses (motif scanning and homology searches).
#'
#' @param src Character. The file path to a single GenBank source file.
#' @param release_dir Character. The root directory where the generated 
#'   \code{builder.gcproj} file will be saved.
#' @param workdir Character. The working directory where extracted FASTA 
#'   files will be written.
#' @param batch_process Logical. If \code{TRUE}, the function extracts the 
#'   model accession ID to create a specific subdirectory inside 
#'   \code{release_dir} for organized batch storage. If \code{FALSE}, files 
#'   are saved directly into the root of \code{release_dir}.
#'
#' @return Returns \code{invisible(NULL)}. Upon execution, it generates the 
#'   following local files:
#'   \itemize{
#'     \item \strong{Project File:} \code{builder.gcproj} containing the 
#'       parsed GenBank project data.
#'     \item \strong{upstream_locis.fasta:} TSS (Transcription Start Site) 
#'       upstream sequences, prepared for downstream motif site scanning.
#'     \item \strong{proteins.fasta:} Genomic protein sequences, prepared 
#'       for enzyme and Transcription Factor (TF) annotation via DIAMOND BLASTp.
#'   }
#'
#' @keywords internal
#' @export
const make_genbank_proj_file = function(src, release_dir, 
                                        workdir = "", 
                                        batch_process = FALSE) {

    let gb_src = load_genbanks(src) |> as.vector();
    let proj = project::new(gb_src);
    let model_id = ifelse(batch_process, model_accession_id(gb_src), ""); 
    let proj_file = file.path(release_dir, model_id, "builder.gcproj");

    # save the ncbi genbank project data as local file
    project::save(proj, file = proj_file);
    # export work files into corresponding model dir
    workdir <- file.path(workdir, model_id);

    # write the TSS upstream and genomics protein fasta sequence
    # as file for the downstream annotation search
    # TSS upstream site for run motif site scanning
    # genomics protein used for make enzyme and TF annotation via the diamond blastp
    # search
    write.fasta(tss_upstream(proj), file = file.path(workdir, "upstream_locis.fasta")); 
    workflow::save_proteins(proj, file = file.path(workdir,"proteins.fasta"));
}