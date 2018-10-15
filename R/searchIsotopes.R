# searchIsotopes
#' Target isotopes search
#'
#' This function uses annotation results of an unlabelled sample to search
#' for labelled compounds in a labelled sample.
#'
#' @param results annotation results for an unlabelled sample. Output of
#' identification functions (i.e. idPOS$results).
#' @param MS1 Data frame with at least three columns: m.z, RT, int. Peak list
#' for the labelled sample. Output of \link{dataProcessing} function
#' (MS$peaklist).
#' @param label isotope employed for the experiment. It can be "13C" or "D".
#' @param adductsTable adducts table employed for lipids annotation.
#' @param rttol rt window in seconds.
#' @param ppm mass error tolerance.
#'
#' @return List with the isotopes for each compound in the results data frame.
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
searchIsotopes <- function(results, MS1, label,
                           adductsTable = LipidMS::adductsTable,
                           rttol = 10, ppm = 10){
  targets <- getInclusionList(results = results, adductsTable = adductsTable)
  comp <- sapply(targets$Formula, CHNOSZ::makeup)
  if("C" %in% names(comp[[1]])){
    targets$C <- unlist(lapply(comp, function(x) x["C"]))
    targets$H <- unlist(lapply(comp, function(x) x["H"]))
  }
  if (label == "13C"){
    massdiff <- 1.0033548
    label <- "C"
  } else if (label == "D"){
    massdiff <- 1.006277
    label <- "H"
  }
  isotopes <- apply(targets, 1, function(x){
    n <- adductsTable[adductsTable[,"adduct"] == as.character(x["Adduct"]), "n"]
    charge <- adductsTable[adductsTable[,"adduct"] == as.character(x["Adduct"]),
                           "charge"]
    mdiff <- adductsTable[adductsTable[,"adduct"] == as.character(x["Adduct"]),
                          "mdiff"]
    mz <- (n*as.numeric(x["Mn"])+as.numeric(mdiff))/abs(charge)
    rt <- as.numeric(x["RT"])
    top <- as.numeric(x[label])
    intensities <- vector()
    mzs <- vector()
    rts <- vector()
    for (i in 0:top){
      if (nrow(MS1[abs(MS1$RT - rt) < rttol,]) > 0){
        subsetMS1 <- MS1[abs(MS1$RT - rt) < rttol,]
        m <- mzMatch(mz+i*massdiff, subsetMS1[abs(subsetMS1$RT - rt) <
                                                         rttol, "m.z"], ppm)
        if (length(m) > 0){
          m <- m[seq(1, length(m), 2)]
          sel <- subsetMS1[m,]
          int <- sel[order(abs(sel$RT - rt)),"int"][1]
          RT <- sel[order(abs(sel$RT - rt)),"RT"][1]
          MZ <- sel[order(abs(sel$RT - rt)),"m.z"][1]
        } else {
          int <- 0
          RT <- rt
          MZ <- mz+i*massdiff
        }
      } else {
        int <- 0
        RT <- rt
        MZ <- mz+i*massdiff
      }
      intensities <- append(intensities, int)
      mzs <- append(mzs, MZ)
      rts <- append(rts, RT)
    }
    return(data.frame(Name = paste(x["Name"], " [M+", 0:top, "]",
                                               sep = ""),
                      Adduct = rep(x["Adduct"], length(mzs)),
                      Formula = rep(x["Formula"], length(mzs)),
                      m.z = mzs, RT = rts,  int = intensities,
                      stringsAsFactors = F))
  })
    return(isotopes)
}

