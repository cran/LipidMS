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
#' @param rawData raw scans data. Output of \link{dataProcessing} function
#' (MS1$rawData).
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied.
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
#' library(LipidMSdata)
#' dbs <- assignDB()
#'
#' candidates <- findCandidates(MS1 =MS1_neg$peaklist,
#' db = dbs$pgdb, ppm = 10, rt = c(0, 2000), adducts = c("M-H"),
#' rttol = 10, rawData = MS1_neg$rawScans, coelCutoff = 0.8)
#'
#'
#' # If any adduct is not in the adductsTable, it can be added:
#'
#' adductsTable2 <- rbind(adductsTable,
#' c(adduct = "M+HCOO", mdiff = 44.9982, n = 1, charge = -1))
#' dbs <- assignDB()
#' dbs$adductsTable <- adductsTable2
#'
#' candidates <- findCandidates(MS1 = MS1_neg$peaklist,
#' db = dbs$pgdb, ppm = 10, rt = c(0, 2000), adducts = c("M-H", "M+HCOO"),
#' rttol = 10, rawData = MS1_neg$rawScans, coelCutoff = 0.8)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
findCandidates <- function(MS1, db, ppm, rt,
                            adducts, rttol = 3, dbs, rawData = data.frame(),
                           coelCutoff = 0){
  adductsTable <- dbs$adductsTable
  count <- 0
  candidates <- vector()
  for (i in 1:length(adducts)){
    ad <- adductsTable[adductsTable$adduct == adducts[i],]
    prec <- findPrecursor(MS1, db, ppm, massdif = ad$mdiff, rt, ad$n, ad$charge)
    if (nrow(prec) > 0){
      prec <- cbind(prec, adducts = as.vector(adducts[i]))
      if (count == 0){
        candidates <- rbind(candidates, prec)
        count <- 1
      } else {
        candidates <- crossAdducts(df1 = candidates, df2 = prec, rttol = rttol,
                                   rawData = rawData, coelCutoff = coelCutoff)
      }
    }
  }
  if(length(candidates) > 0){
    candidates <- filtrateAdducts(df = candidates)
  }
  rownames(candidates) <- c()
  if (is.vector(candidates)){
    return(data.frame())
  } else {
    return(candidates)
  }
}

