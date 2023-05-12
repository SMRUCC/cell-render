import {jQuery, Html, http} from "webKit";

setwd(@dir);

var groups = JSON.parse(readText('./taxonomics_group.json'));
var local = http.cache("./cache/");

// console.log(groups);

for(var group in groups) {
    var html = jQuery.load(sprintf('https://regprecise.lbl.gov/%s', group.id), proxy = local);
    var tbl = html[".stattbl"]
    var body = tbl["tbody"]
    var rows = body["tr"]
    var taxonomics = lapply(rows, function(r) {
        var cells = r["td"]
        var name = cells[1]
        
        name = name.innerHTML;
        
        var id = Html.link(name);
        
        name = Html.plainText(name);
        
        return {
            id: id, 
            name: name
        }
    });

    console.log(group);
    console.log(taxonomics);
    
    throw 'just stop for debug'
}