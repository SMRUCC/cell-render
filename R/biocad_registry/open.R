#' open mysql connection to the cad registry
#' 
const open_registry = function(user,passwd, host = "localhost",port = 3306) {
    require(graphQL);

    # load mysql driver
    imports "mysql" from "graphR";
    # create registry connection
    mysql::open(
        user_name = user,
        password = passwd,
        dbname = "cad_registry",
        host  = host,
        port = port
    );
}