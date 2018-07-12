# findCandidates
#' Search of lipid candidates of a certain class.
#'
#' Search of lipid candidates from a peaklist based on a set of expected
#' adducts.
#'
#' @param MS1 peaklist of the MS function. Data frame with 3 columns: m.z, RT
#' (in seconds) and int (intensity).
#' @param db database (i.e. pcdb, dgdb, etc.). Data frame with at least 2
#' columns: Mass (exact mass) and total (total number of carbons and double
#' bound of the FA chains, i.e. "34:1").
#' @param ppm m/z tolerance in ppm.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts character vector containing the expected adducts to search
#' for (i.e. "M+H", "M+Na", "M-H", etc.). See details.
#' @param rttol rt tolerance in seconds to match adducts.
#' @param dbs list of data bases required for the annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be changed. If data bases have
#' been customized using \link{createLipidDB}, they also have to be modified
#' here.
#'
#' @return Data frame with the found candidates. It contains 6 columns: m.z,
#' RT, int (from the peaklist data.frame), ppms, cb (total number of carbons and
#' double bounds of the FA chains) and adducts.
#'
#' @details \link{findCandidates} looks for matches between the m/z of the
#' MS1 peaklist and the expected m/z of the candidates in the database for each
#' adduct. If several adducts are expected, results are combined.
#'
#' Adducts allowed are contained in adductsTable data frame, which can be
#' modified if required (see \link{adductsTable}).
#'
#' @examples
#' \donttest{
#' dbs <- list(cerdb = LipidMS::cerdb, adductsTable = LipidMS::adductsTable)
#'
#' findCandidates(MS1 = LipidMS::mix_neg_fullMS,
#' db = dbs[["cerdb"]], ppm = 5, rt = c(min(MS1$RT), max(MS1$RT)),
#' adducts = c("M-H", "M-H-H2O", "M+CH3COO"), rttol = 3, dbs)
#'
#' #If any adduct is not in the adductsTable, it can be added:
#'
#' adductsTable2 <- rbind(LipidMS::adductsTable,
#' c(adduct = "M+HCOO", mdiff = 44.99820, n = 1, charge = -1))
#' dbs <- list(cerdb = LipidMS::cerdb, adductsTable = adductsTable2)
#'
#' findCandidates(MS1 = LipidMS::mix_neg_fullMS,
#' db = dbs[["cerdb"]], ppm = 5, rt = c(min(MS1$RT), max(MS1$RT)),
#' adducts = c("M-H", "M-H-H2O", "M+HCOO"), rttol = 3, dbs)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
findCandidates <- function(MS1, db, ppm = 10, rt = c(min(MS1$RT), max(MS1$RT)),
                           adducts, rttol = 3, dbs){
  adductsTable <- dbs[["adductsTable"]]
  count <- 0
  candidates <- vector()
  for (i in 1:length(adducts)){
    ad <- adductsTable[adductsTable$adduct == adducts[i],]
    prec <- findPrecursor(MS1, db, ppm, ad$mdiff, rt, ad$n, ad$charge)
    if (nrow(prec) > 0){
      if (count == 0){
        candidates <- rbind(candidates, cbind(prec, adducts = adducts[i]))
        count <- 1
      } else {
        candidates <- joinAdducts(candidates, prec,
                                  rttol, "", adducts[i])
      }
    }
  }
  if (is.vector(candidates)){
    return(data.frame())
  } else {
    return(candidates)
  }
}
