#' Register the Virtual Cell Network Modelling Workflows
#'
#' Defines and registers the sequential annotation and modelling workflow
#' steps that constitute the virtual cell graph construction pipeline.
#' Each step is hooked into the \pkg{WorkflowRender} engine and will be
#' executed in the order they are registered.
#'
#' @details
#' The complete workflow sequence is:
#' \enumerate{
#'   \item \strong{make_genbank_proj} \cr
#'     Extract protein sequences, TSS upstream regions, genomic context data,
#'     and taxonomy information from the GenBank source file(s).
#'   \item \strong{tfbs_motif_scanning} \cr
#'     Run Transcription Factor Binding Site (TFBS) motif scanning on the
#'     TSS upstream region sequences extracted in step 1.
#'   \item \strong{make_diamond_hits} \cr
#'     Run DIAMOND BLASTP homology search for the protein sequences against
#'     EC number, subcellular location, and transcription factor reference
#'     databases.
#'   \item \strong{make_terms} \cr
#'     Score the BLASTP results and assign protein ontology terms (EC numbers,
#'     subcellular localizations, transcription factor annotations).
#'   \item \strong{make_TRN} \cr
#'     Construct the Transcription Regulation Network (TRN) by combining the
#'     TFBS motif scan results with the transcription factor BLASTP hits.
#'   \item \strong{build_project} \cr
#'     Compile the TRN and metabolic network annotations into a final virtual
#'     cell model and export as GCMarkup XML.
#' }
#'
#' @return \code{invisible(NULL)}. This function is called for its side effect
#'   of registering workflow steps with the \pkg{WorkflowRender} engine.
#'
#' @seealso \code{\link{make_genbank_proj}}, \code{\link{tfbs_motif_scanning}},
#'   \code{\link{make_diamond_hits}}, \code{\link{make_terms}},
#'   \code{\link{make_TRN}}, \code{\link{build_project}}
#'
#' @examples
#' \dontrun{
#' # Register all workflow steps
#' annotation_workflow()
#'
#' # Then run the workflow
#' WorkflowRender::run(registry = CellRender::annotation_workflow)
#' }
#'
#' @app annotation_workflow
#' @export
const annotation_workflow = function(debug) {
    print("Config the cellular graph network annotation & modelling workflow...");

    # [@app "xxx"] defines the workflow module, and "xxx" is the workflow module name
    # current workflow sequence is:
    #
    # step1. make_genbank_proj - extract protein sequence and TSS upstream region, genomics context data, taxonomy information from the genbank file
    # step2. tfbs_motif_scanning - run TFBS motif site scanning on the TSS upstream region
    # step3. make_diamond_hits - run diamond blastp search for the protein sequence
    # step5. make_terms - run blastp result scoreing and assigned protein ontology term
    # step6. make_TRN - create transcription regulation network based on the TFBS site scaned from the TSS upstream region and TF blastp result 
    # step7. build_project - compile the TRN and metabolic network as virtual cell model

    # steps for make genbank annotation project
    WorkflowRender::hook(CellRender::make_genbank_proj);
    WorkflowRender::hook(CellRender::tfbs_motif_scanning);
    # blastp search and then make annotation terms
    WorkflowRender::hook(CellRender::make_diamond_hits);
    WorkflowRender::hook(CellRender::make_terms);
    WorkflowRender::hook(CellRender::make_TRN);

    # steps for create virtual cell model from genbank annotation project file
    WorkflowRender::hook(CellRender::build_project);

    WorkflowRender::summary();

    if (length(debug) > 0) {
        WorkflowRender::definePipeline(debug);
    }    

    invisible(NULL);
} 
