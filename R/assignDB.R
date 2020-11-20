# assignDB
#' load LipidMS default data bases
#'
#' load all LipidMS default data bases required to run identification functions.
#'
#' @return list of data frames
#'
#' @examples
#' \donttest{
#' dbs <- assignDB()
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
assignDB <- function(){
  dbs <- list()
  dbs[["cerdb"]] <- LipidMS::cerdb
  dbs[["smdb"]] <- LipidMS::smdb
  dbs[["fadb"]] <- LipidMS::fadb
  dbs[["hfadb"]] <- LipidMS::hfadb
  dbs[["carnitinedb"]] <- LipidMS::carnitinesdb
  dbs[["lysopadb"]] <- LipidMS::lysopadb
  dbs[["lysopedb"]] <- LipidMS::lysopedb
  dbs[["lysopgdb"]] <- LipidMS::lysopgdb
  dbs[["lysopidb"]] <- LipidMS::lysopidb
  dbs[["lysopsdb"]] <- LipidMS::lysopsdb
  dbs[["lysopcdb"]] <- LipidMS::lysopcdb
  dbs[["mgdb"]] <- LipidMS::mgdb
  dbs[["CEdb"]] <- LipidMS::CEdb
  dbs[["sphdb"]] <- LipidMS::sphdb
  dbs[["sphPdb"]] <- LipidMS::sphPdb
  dbs[["fahfadb"]] <- LipidMS::fahfadb
  dbs[["pcdb"]] <- LipidMS::pcdb
  dbs[["pedb"]] <- LipidMS::pedb
 dbs[["pgdb"]] <- LipidMS::pgdb
  dbs[["pidb"]] <- LipidMS::pidb
  dbs[["psdb"]] <- LipidMS::psdb
  dbs[["padb"]] <- LipidMS::padb
  dbs[["dgdb"]] <- LipidMS::dgdb
  dbs[["tgdb"]] <- LipidMS::tgdb
  dbs[["cldb"]] <- LipidMS::cldb
  dbs[["badb"]] <- LipidMS::badb
  dbs[["baconjdb"]] <- LipidMS::baconjdb
  dbs[["nlsphdb"]] <- LipidMS::nlsphdb
  dbs[["adductsTable"]] <- LipidMS::adductsTable
  return(dbs)
}
