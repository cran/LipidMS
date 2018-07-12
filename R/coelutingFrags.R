# coelutingFrags
#' Coeluting fragments extraction
#'
#' Given a RT and a list of peaks, this function subsets all coeluting fragments
#' within a rt windows. It is used by identification functions to extract
#' coeluting fragments from high energy functions for candidate precursor ions.
#'
#' @param rt numeric vector indicating candidates RT. This information comes
#' from the output of \link{findCandidates} (candidates$RT).
#' @param df data frame containing the peaks to subset (MSMS).
#' @param rttol rt window in seconds.
#'
#' @return List of data frames with the coeluting fragments for each candidate.
#'
#' @examples
#' \donttest{
#' dbs <- list(pgdb = LipidMS::pgdb, lysopgdb = LipidMS::lysopgdb,
#' fadb = LipidMS::fadb, adductsTable = LipidMS::adductsTable)
#'
#' candidates <- findCandidates(MS1 = LipidMS::mix_neg_fullMS, dbs[["pgdb"]],
#' ppm = 10, rt = c(min(MS1$RT), max(MS1$RT)), adducts = c("M-H"),
#' rttol = 3, dbs)
#'
#' MSMS <- rbind(LipidMS::mix_neg_Ce20, LipidMS::mix_neg_Ce40)
#' coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
coelutingFrags <- function(rt, df, rttol){
  psubset <- lapply(rt, function(x){
    peaks <- which(df[,"RT"] <= x+(rttol/2) & df[,"RT"] >= x-(rttol/2))
    df2 <- df[peaks, c("m.z", "int", "RT")]
    return(df2)
  })
  return(psubset)
}
