imports "umap" from "MLkit";

setwd(@dir);

let expr = read.csv("./metabolites.csv", row.names = 1, check.names = FALSE);
let result = umap(expr, dimension = 3, numberOfNeighbors = 128);
let manifold = as.data.frame(result$umap, labels = result$labels);

write.csv(manifold, file = "./umap3.csv");