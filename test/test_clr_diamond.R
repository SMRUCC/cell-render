require(CellRender);

let q = read.fasta( "F:\ecoli\29.0\data\protseq.fsa");
let s = read.fasta("G:\cell-render\data\KO.fasta");
let result = clr_diamond (q,s);

write.csv(result , file = here("ecoli_ko.csv"));
