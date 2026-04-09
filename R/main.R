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
#' @param src \code{character}. The file path to the NCBI GenBank assembly 
#'   source file.
#' @param outputdir \code{character}. The directory path for storing result 
#'   files and temporary workspace. Defaults to \code{NULL}, which will use 
#'   the parent directory of the \code{src} file.
#' @param name \code{character}. An optional custom name for the virtual cell 
#'   model. If \code{NULL}, a default name may be assigned by the rendering engine.
#' @param up_len \code{integer}. The length (in base pairs) of the upstream 
#'   region from the Transcription Start Site (TSS) to be extracted. 
#'   Defaults to 150.
#' @param localdb \code{character}. The file path to a local annotation 
#'   database. Defaults to \code{NULL}, which falls back to the internal 
#'   package data directory.
#' @param diamond \code{character}. The file path to the DIAMOND executable 
#'   (used for fast sequence alignment). Defaults to the system path found 
#'   via \code{Sys.which("diamond")}.
#' @param domain \code{character}. The biological domain of the input organism. 
#'   Must be one of \code{"bacteria"}, \code{"plant"}, \code{"animal"}, or 
#'   \code{"fungi"}. Note: If multiple values are provided, only the first 
#'   element will be used.
#' @param builds \code{character}. A character vector specifying which 
#'   annotation modules to build. Supported options are \code{"TRN_network"} 
#'   (Transcriptional Regulatory Network) and \code{"Metabolic_network"}.
#' @param n_threads \code{integer}. The number of CPU threads to allocate 
#'   for parallel processing. Defaults to 32.
#'
#' @return This function returns \code{invisible(NULL)}. However, as a side 
#'   effect, it generates a \code{release} subdirectory within the 
#'   \code{outputdir} containing:
#'   \itemize{
#'     \item \code{builder.gcproj}: The project configuration file.
#'     \item \code{model.xml}: The core virtual cell graph model file.
#'   }
#'   Additionally, an HTML or PDF report summarizing the virtual cell 
#'   modeling results is typically generated in the output directory by 
#'   the workflow renderer.
#'
#' @examples
#' \dontrun{
#' # Run the workflow for a bacterial genome
#' modelling_cellgraph(
#'   src = "path/to/bacteria_genome.gb",
#'   outputdir = "path/to/output_workspace",
#'   name = "E_coli_Model",
#'   domain = "bacteria",
#'   builds = c("TRN_network", "Metabolic_network"),
#'   n_threads = 16
#' )
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
                                     n_threads = 32, 
                                     debug = c()) {

    let batch_process as boolean = dir.exists(src); 

    WorkflowRender::init_context(outputdir || dirname(src));
    WorkflowRender::set_config(list(
        src        = normalizePath(src),
        localdb    = localdb || normalizePath(@datadir),
        up_len     = up_len,
        diamond    = unlist(diamond),
        n_threads  = n_threads,
        domain     = .Internal::first(domain),
        builds     = builds,
        release    = file.path(workdir_root(), "release"),
        gem_libout = file.path(workdir_root(), "GEMs"),
        proj_file  = file.path(workdir_root(), "release", "builder.gcproj"),
        model_file = file.path(workdir_root(), "release", "model.xml"),
        vcell_name = name,
        # the input source filesystem handle is a directory
        # that contains multiple genbank assembly files
        # run this workflow in batch mode 
        batch_process = batch_process,
        enable_blastp_cache = enable_blastp_cache,
        gems_library_mode = gems_library_mode,
        enzyme_fuzzy = enzyme_fuzzy
    ));

    if (length(debug) > 0) {
        WorkflowRender::definePipeline(debug);
    }    

    WorkflowRender::summary();
    WorkflowRender::run(registry = CellRender::annotation_workflow);
    WorkflowRender::finalize();

    invisible(NULL);
}

