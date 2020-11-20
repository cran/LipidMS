# coelutingFrags
#' Coeluting fragments extraction
#'
#' Given a RT and a list of peaks, this function subsets all coeluting fragments
#' within a rt windows. It is used by identification functions to extract
#' coeluting fragments from high energy functions for candidate precursor ions.
#'
#' @param precursors candidates data frame. Output of \link{findCandidates}.
#' @param products peaklist for MS2 function (MSMS).
#' @param rttol rt window in seconds.
#' @param rawData raw scans data. Output of \link{dataProcessing} function
#' (MSMS$rawData).
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied.
#'
#' @return List of data frames with the coeluting fragments for each candidate.
#'
#' @examples
#' \donttest{
#' library(LipidMSdata)
#' dbs <- assignDB()
#'
#' candidates <- findCandidates(MS1 = MS1_neg$peaklist,
#' db = dbs$pgdb, ppm = 10, rt = c(0, 2000), adducts = c("M-H"),
#' rttol = 10, dbs = dbs, rawData = MS1_neg$rawScans, coelCutoff = 0.8)
#'
#' MSMS <- rbind(MSMS1_neg$peaklist, MSMS2_neg$peaklist)
#' rawData <- rbind(MS1_neg$rawScans, MSMS1_neg$rawScans, MSMS2_neg$rawScans)
#' coelfrags <- coelutingFrags(candidates, MSMS, rttol = 10, rawData = rawData,
#' coelCutoff = 0.8)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
coelutingFrags <- function(precursors, products, rttol, rawData = data.frame(),
                           coelCutoff = 0){
  psubset <- apply(precursors, 1, function(x){
    peaks <- which(products[,"RT"] <= as.numeric(x["RT"])+(rttol/2) & products[,"RT"] >=
                     as.numeric(x["RT"])-(rttol/2))
    if(length(peaks) > 0){
      df <- products[peaks,]
      scores <- coelutionScore(as.character(x["peakID"]), df$peakID, rawData)
      df$coelScore <- scores
      df <- df[scores >= coelCutoff,]
      return(df)
    } else {
      return(data.frame())
    }
  })
  coelfrags <- lapply(psubset, function(x) {
    if (nrow(x) > 0){
      x <- x[!is.na(x$m.z),]
    } else {data.frame()}
  })
  return(coelfrags)
}


