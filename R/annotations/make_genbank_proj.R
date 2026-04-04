imports "project" from "CellRender";
imports "workflow" from "CellRender";

#' Initialize a GenBank Project and Export Downstream Sequence Files
#'
#' @description
#' Loads GenBank source files specified in the application configuration,
#' initializes a new `project` object via `CellRender`, and saves the 
#' resulting project state to disk. Additionally, it extracts and exports 
#' essential FASTA files required for downstream bioinformatics analyses,
#' such as motif scanning and homology searches.
#'
#' @param app A CellRender workflow app object. Primarily used in this 
#'   function to resolve working file paths via \code{workfile()}.
#' @param context A CellRender workflow context object. 
#'   \emph{Note:} Currently unused in this specific function but required 
#'   by the standard workflow execution signature.
#'
#' @return 
#' This function is called for its side effects and returns \code{NULL} 
#' invisibly. Upon execution, it generates the following local files:
#' \itemize{
#'   \item \strong{Project File:} The serialized project object (path defined 
#'     by \code{get_config("proj_file")}).
#'   \item \strong{upstream_locis.fasta:} Contains TSS (Transcription Start Site) 
#'     upstream sequences, prepared for downstream motif site scanning.
#'   \item \strong{proteins.fasta:} Contains genomic protein sequences, prepared 
#'     for enzyme and Transcription Factor (TF) annotation via Diamond BLASTp.
#' }
#'
#' @details
#' The function executes a sequential pipeline:
#' \enumerate{
#'   \item \strong{Data Ingestion:} Reads GenBank files from the path 
#'     specified in \code{get_config("src")} using \code{load_genbanks()}.
#'   \item \strong{Object Creation:} Wraps the raw GenBank data into a 
#'     standardized \code{project} object.
#'   \item \strong{State Persistence:} Saves the core project object locally 
#'     to avoid redundant parsing in future workflow steps.
#'   \item \strong{Sequence Export:** Formats and writes specific sequence 
#'     subsets (upstream loci and proteins) into standardized FASTA files 
#'     tailored for specific downstream tools.
#' }
#'
#' @importFrom CellRender project
#' @importFrom CellRender workflow
#'
#' @examples
#' # Typically called internally by the CellRender workflow engine:
#' # make_genbank_proj(current_app, current_context)
#'
[@app "make_genbank_proj"]
const make_genbank_proj = function(app, context) {
    let batch_process = as.logical(get_config("batch_process"));
    let release_dir = get_config("release");
    let workdir = WorkflowRender::workspace(app);

    if (batch_process) {  
        let genbank_files = list.files(get_config("src"), 
                                pattern = c("*.gb","*.gbk","*.gbff"));

        message(`Build virtualcell community model based ${length(genbank_files)} genbank source files!`);

        for(let file in tqdm(genbank_files)) {
            file |> make_genbank_proj_file(
                release_dir = release_dir,
                workdir = workdir,
                batch_process = TRUE
            );
        }
    } else {
        # input source is a single genbank file
        let gb_src = get_config("src");        

        gb_src |> make_genbank_proj_file(
            release_dir = release_dir,
            workdir = workdir,
            batch_process = FALSE
        );
    }

    invisible(NULL);
}

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