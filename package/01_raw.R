imports "visualPlot" from "visualkit";
imports "clustering" from "MLkit";
imports "charts" from "R.plot";

let hist_samples as function(workspace, raw) {
    let outputdir  = `${workspace$dirs$stat}`;
    let sampleinfo = sampleclass(as.data.frame(workspace$sample_info));
    let d = raw 
    :> t 
    :> dist 
    :> hclust
    ;

    print(d);
    str(sampleinfo);

    d :> plot(
        class       = sampleinfo, 
        size        = [2700, 4000], 
        padding     = "padding: 200px 400px 200px 200px;", 
        axis.format = "G3",
        links       = "stroke: darkblue; stroke-width: 8px; stroke-dash: dash;",
        pt.color    = "gray",
        label       = "font-style: normal; font-size: 13; font-family: Bookman Old Style;",
        ticks       = "font-style: normal; font-size: 10; font-family: Bookman Old Style;"
    )
    :> save.graphics(`${outputdir}/hclust_samples.png`)
    ;
}

let sampleclass as function(sampleinfo) {
    let ID     = sampleinfo[, "ID"];
    let colors = "#" & sampleinfo[, "color"];

    print("previews of your raw sample information:");
    print(head(sampleinfo));

    lapply(1:length(ID), i -> colors[i], names = i -> ID[i]);
}