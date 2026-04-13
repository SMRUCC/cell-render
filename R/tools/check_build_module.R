
const check_build_module = function(flag) {
    let builds = get_config("builds");
    let check = any(tolower(flag) == tolower(builds));

    return(check);
}