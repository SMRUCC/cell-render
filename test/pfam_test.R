require(CellRender);
require(GCModeller);

imports "annotation.workflow" from "seqtoolkit";
imports "proteinKit" from "seqtoolkit";

# CellRender::pfam_diamond(
#     proteins = "./models/tmp/workflow_tmp/make_genbank_proj/NZ_DS499581/proteins.fasta", 
#     workdir = "./models/tmp/workflow_tmp/make_genbank_proj/NZ_DS499581", 
#     diamond = Sys.which("diamond"));

let pfam = "\\192.168.3.15\sda\2026\Tryptophan_20260520\models\tmp\workflow_tmp\make_genbank_proj\NZ_DS499581\proteins.tsv";

pfam = read_m8(pfam);
pfam = proteinKit::analysis_domains(pfam);

write.csv(pfam, file = "Z:/test.csv");