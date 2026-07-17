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
#'   \code{release_dir} for organizing multiple models. If \code{FALSE},
#'   files are written directly to the provided directories without
#'   subdirectory nesting.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects:
#'   \itemize{
#'     \item Creates a \code{builder.gcproj} project file in the release directory
#'       (under a model-ID subdirectory in batch mode).
#'     \item Exports \code{upstream_locis.fasta} containing TSS upstream region
#'       sequences for downstream TFBS motif scanning.
#'     \item Exports \code{proteins.fasta} containing protein sequences for
#'       downstream DIAMOND BLASTP annotation searches.
#'   }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Reads the GenBank source file and creates a \code{project} object.
#'   \item Extracts the model accession ID (in batch mode) for directory naming.
#'   \item Saves the project file to the release directory.
#'   \item Extracts and writes TSS upstream sequences (\code{upstream_locis.fasta})
#'     for transcription factor binding site motif scanning.
#'   \item Extracts and writes protein sequences (\code{proteins.fasta}) for
#'     DIAMOND BLASTP homology annotation.
#' }
#'
#' @seealso \code{\link{make_genbank_proj}} which calls this function,
#'   \code{\link{model_accession_id}} for accession ID extraction.
#'
#' @examples
#' \dontrun{
#' # Process a single GenBank file
#' make_genbank_proj_file(
#'   src = "GCF_000123456.1.gbff",
#'   release_dir = "release",
#'   workdir = "workspace/make_genbank_proj",
#'   batch_process = FALSE
#' )
#'
#' # Process in batch mode (creates subdirectory by accession ID)
#' make_genbank_proj_file(
#'   src = "GCF_000123456.1.gbff",
#'   release_dir = "release",
#'   workdir = "workspace/make_genbank_proj",
#'   batch_process = TRUE
#' )
#' }
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

    return(dirname(proj_file));
}