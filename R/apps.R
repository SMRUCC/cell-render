#' Register the virtual cell network modelling workflows
#' 
#' @return this function has no return value
#' 
const annotation_workflow = function() {
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

    invisible(NULL);
} 
