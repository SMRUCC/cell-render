#' Run the Virtual Cell Graph Modeling Workflow
#' 
#' Executes the end-to-end workflow for constructing a virtual cell graph model
#' from an NCBI GenBank assembly source file. This function initializes the 
#' environment, configures parameters, runs the annotation workflow, and 
#' generates the corresponding project files and reports.
#' 
#' @details 
#' This function acts as a high-level wrapper around the \code{WorkflowRender} 
#' and \code{CellRender} engines. It sets up a temporary workspace (if not 
#' provided), resolves dependencies (like the DIAMOND executable), and triggers 
#' the specific annotation modules requested by the user.
#'
#' The complete workflow pipeline includes:
#' \enumerate{
#'   \item \strong{make_genbank_proj}: Extract protein sequences, TSS upstream
#'     regions, genomic context, and taxonomy from GenBank source file(s).
#'   \item \strong{tfbs_motif_scanning}: Scan for transcription factor binding
#'     site motifs in upstream regions (conditional on TRN_network module).
#'   \item \strong{make_diamond_hits}: Run DIAMOND BLASTP against EC number,
#'     subcellular location, and transcription factor reference databases.
#'   \item \strong{make_terms}: Parse BLASTP results and assign annotation terms
#'     for metabolic, transmembrane, and transcription regulation networks.
#'   \item \strong{make_TRN}: Build the transcription regulation network from
#'     TFBS scan results and TF annotations (conditional on TRN_network module).
#'   \item \strong{build_project}: Compile all networks into a virtual cell
#'     model and export as GCMarkup XML.
#' }
#'
#' @param src \code{character}. The file path to the NCBI GenBank assembly 
#'   source file. In batch mode, this can be a directory containing multiple
#'   GenBank files.
#' @param outputdir \code{character}. The directory path for storing 
#'   result files. Defaults to a temporary directory created by
#'   \code{tempdir()} if not specified.
#' @param name \code{character}. The name to assign to the virtual cell model.
#'   Only used in single-project mode. Defaults to \code{"virtual_cell"}.
#' @param batch_process \code{logical}. If \code{TRUE}, treats \code{src} as
#'   a directory and processes all GenBank files within it. If \code{FALSE},
#'   processes a single GenBank file. Defaults to \code{FALSE}.
#' @param debug \code{character vector}. Optional vector of workflow step names
#'   to enable for debugging. When provided, only the specified steps will be
#'   executed via \code{WorkflowRender::definePipeline()}. Defaults to an
#'   empty vector (all steps run).
#' @param enable_blastp_cache \code{logical}. If \code{TRUE}, enables caching
#'   of BLASTP search results to avoid redundant computations on re-runs.
#'   Defaults to \code{FALSE}.
#' @param gems_library_mode \code{logical}. If \code{TRUE}, runs in GEMs
#'   library mode, which adjusts the model compilation for genome-scale
#'   metabolic model library generation. Defaults to \code{FALSE}.
#' @param enzyme_fuzzy \code{logical}. If \code{TRUE}, enables fuzzy matching
#'   for enzyme annotation, allowing more permissive EC number assignments.
#'   Defaults to \code{FALSE}.
#'
#' @return \code{invisible(NULL)}. This function is called for its side effects
#'   of creating project files, annotation data, and compiled virtual cell
#'   model XML files in the output directory.
#'
#' @seealso \code{\link{annotation_workflow}} for the workflow registration
#'   function, \code{\link{make_genbank_proj}}, \code{\link{make_diamond_hits}},
#'   \code{\link{make_terms}}, \code{\link{make_TRN}}, \code{\link{build_project}}
#'   for individual workflow steps.
#'
#' @examples
#' \dontrun{
#' # Single GenBank file mode
#' run(src = "data/GCF_000005845.2_ASM584v2_genomic.gbff",
#'     name = "E_coli_K12",
#'     outputdir = "results/E_coli")
#'
#' # Batch mode for multiple GenBank files
#' run(src = "data/genbank_assemblies/",
#'     batch_process = TRUE,
#'     outputdir = "results/batch")
#' }
#'
#' @export
const modelling_cellgraph = function(src, outputdir = NULL, 
                                     name = NULL,
                                     up_len = 150, 
                                     localdb = NULL, 
                                     diamond = Sys.which("diamond"), 
                                     domain = c("bacteria", "plant", "animal", "fungi"),
                                     builds = c("TRN_network","Metabolic_network"),
                                     enable_blastp_cache = FALSE,
                                     enzyme_fuzzy = FALSE,
                                     gems_library_mode = TRUE,
                                     batch_mode = c("batch","sequential"),
                                     n_threads = 32, 
                                     debug = c()) {

    let batch_process as boolean = dir.exists(src); 
    let workdir_root = outputdir || dirname(src);
    let args = list(
        src        = normalizePath(src),
        localdb    = localdb || normalizePath(@datadir),
        up_len     = up_len,
        diamond    = unlist(diamond),
        n_threads  = n_threads,
        domain     = .Internal::first(domain),
        builds     = builds,
        release    = file.path(workdir_root, "release"),
        gem_libout = file.path(workdir_root, "GEMs"),
        proj_file  = file.path(workdir_root, "release", "builder.gcproj"),
        model_file = file.path(workdir_root, "release", "model.xml"),
        vcell_name = name,
        # the input source filesystem handle is a directory
        # that contains multiple genbank assembly files
        # run this workflow in batch mode 
        batch_process = batch_process,
        enable_blastp_cache = enable_blastp_cache,
        gems_library_mode = gems_library_mode,
        enzyme_fuzzy = enzyme_fuzzy,
        batch_mode = batch_mode
    );

    WorkflowRender::init_context(workdir_root );
    WorkflowRender::set_config(args);

    if (batch_mode == "sequential") {
        sequential_batch(src, 
            outputdir = workdir_root , 
            args = args);
    } else {
        if (length(debug) > 0) {
            WorkflowRender::definePipeline(debug);
        }    

        # config run workflow
        WorkflowRender::run(registry = CellRender::annotation_workflow);
        WorkflowRender::finalize();
    }

    invisible(NULL);
}

