imports "http" from "webKit";

const taxonomic_group.create = function(name, note = "") {
    const base = getOption("biocad");
    const url  = `${base}/registry/put/taxonomic/`;
    const resp = requests.post(url, list(name, note)) |> http::content();

    str(resp);
}