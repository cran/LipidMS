# crossTables
#' Cross the original MS1 peaklist with the annotation results
#'
#' Cross the original MS1 peaklist with the annotation results.
#'
#' @param MS1 data frame cointaining all peaks from the full MS function. It
#' must have three columns: m.z, RT (in seconds) and int (intensity).
#' @param results data frame. Output of identification functions.
#' @param ppm mass tolerance in ppm.
#' @param rttol rt tolerance to match peaks in seconds.
#' @param adductsTable data frame with the adducts to employ.
#'
#' @return Data frame with 6 columns: m.z, RT, int, LipidMS_id, adduct and
#' confidence level for the annotation. When multiple IDs are proposed for the
#' same feature, they are sorted based on the annotation level.
#'
#' @examples
#' \donttest{results <- idNEG(LipidMS::serum_neg_fullMS, LipidMS::serum_neg_Ce20,
#' LipidMS::serum_neg_Ce40)
#' crossTables(MS1 = LipidMS::serum_neg_fullMS, results = results$results,
#' ppm = 10, rttol = 5, adductsTable = LipidMS::adductsTable)}
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
crossTables <- function(MS1, results, ppm=5, rttol=5,
                        adductsTable = LipidMS::adductsTable){
  confLevels <- LipidMS::confLevels
  rawPeaks <- data.frame(MS1, LipidMS_id = "", Adduct = "",
                         confidenceLevel = "", stringsAsFactors = F)
  Form_Mn <- apply(results, 1, getFormula)
  if (class(Form_Mn) == "matrix"){
    new <- list()
    for (i in 1:ncol(Form_Mn)){
      new[[i]] <- c(Form_Mn[1,i], Form_Mn[2,i])
    }
    Form_Mn <- new
  }
  na <- which(unlist(lapply(Form_Mn, length)) == 0)
  if (length(na) > 0){
    Form_Mn[[na]] <- c(NA, NA)
  }
  Mn <- as.numeric(unlist(lapply(Form_Mn, "[[", 2)))
  for (r in 1:nrow(results)){
    adducts <- unlist(strsplit(as.character(results$Adducts)[r], ";"))
    mzs <- vector()
    for (a in 1:length(adducts)){
      adinfo <- adductsTable[adductsTable$adduct == adducts[a],]
      m.z <- (adinfo$n*Mn[r]+adinfo$mdif)/abs(adinfo$charge)
      mzs <- append(mzs, m.z)
    }
    if (length(mzs) > 0){
      for (m in 1:length(mzs)){
        if (sum(abs(rawPeaks$RT- as.numeric(results$RT[r])) < rttol) > 0){
          rows <- which(abs(rawPeaks$RT - as.numeric(results$RT[r])) < rttol)
          matches <- as.numeric(unlist(sapply(mzs[m], mzMatch,
                                              rawPeaks$m.z[rows], ppm)))
          matches <- matches[seq(1, length(matches), 2)]
        } else {
          matches <- vector()
        }
        if (length(matches) > 0){
          for (i in 1:(length(matches))){
            if (rawPeaks$LipidMS_id[rows[matches[i]]] == ""){
              rawPeaks$LipidMS_id[rows[matches[i]]] <-
                as.character(results$ID[r])
              rawPeaks$confidenceLevel[rows[matches[i]]] <-
                as.character(results$confidenceLevel[r])
              rawPeaks$Adduct[rows[matches[i]]] <-
                as.character(adducts[m])
            } else if (rawPeaks$LipidMS_id[rows[matches[i]]] != "" &&
                        !grepl(results$ID[r],
                               rawPeaks$LipidMS_id[rows[matches[i]]], fixed = T)){
              conf <- confLevels[confLevels$level ==
                                   results$confidenceLevel[r], "order"]
              conf2 <- unlist(strsplit(
                rawPeaks$confidenceLevel[rows[matches[i]]], "\\|"))
              conf2 <- sapply(conf2, function(x)
                confLevels[confLevels$level == x,"order"])
              newposition <- which(conf2 < conf)
              if (length(newposition) == 0){
                newposition <- 0
              }
              if (newposition == 0){
                rawPeaks$LipidMS_id[rows[matches[i]]] <-
                  paste(rawPeaks$LipidMS_id[rows[matches[i]]],
                        as.character(results$ID[r]), sep="|")
                rawPeaks$confidenceLevel[rows[matches[i]]] <-
                  paste(rawPeaks$confidenceLevel[rows[matches[i]]],
                        as.character(results$confidenceLevel[r]), sep="|")
                rawPeaks$Adduct[rows[matches[i]]] <-
                  paste(as.character(rawPeaks$Adduct[rows[matches[i]]]),
                        as.character(adducts[m]), sep="|")
              } else if (newposition == 1){
                rawPeaks$LipidMS_id[rows[matches[i]]] <-
                  paste(as.character(results$ID[r]),
                        rawPeaks$LipidMS_id[rows[matches[i]]], sep="|")
                rawPeaks$confidenceLevel[rows[matches[i]]] <-
                  paste(as.character(results$confidenceLevel[r]),
                        rawPeaks$confidenceLevel[rows[matches[i]]], sep="|")
                rawPeaks$Adduct[rows[matches[i]]] <-
                  paste(as.character(adducts[m]),
                        as.character(rawPeaks$Adduct[rows[matches[i]]]), sep="|")
              } else {
                ids <- unlist(strsplit(as.character(rawPeaks$LipidMS_id[rows[matches[i]]]),
                                       "\\|"))
                conflev <- unlist(strsplit(as.character(rawPeaks$confidenceLevel[rows[matches[i]]]),
                                           "\\|"))
                adduc <- unlist(strsplit(as.character(rawPeaks$Adduct[rows[matches[i]]]),
                                         "\\|"))
                rawPeaks$LipidMS_id[rows[matches[i]]] <-
                  paste(paste(ids[1:(newposition-1)], collapse="|"),
                        as.character(results$ID[r]),
                        paste(ids[newposition:length(ids)], collapse = "|"),
                        sep="|")
                rawPeaks$confidenceLevel[rows[matches[i]]] <-
                  paste(paste(conflev[1:(newposition-1)], collapse="|"),
                        as.character(results$confidenceLevel[r]),
                        paste(conflev[newposition:length(conflev)],
                              collapse = "|"),
                        sep="|")
                rawPeaks$Adduct[rows[matches[i]]] <-
                  paste(paste(adduc[1:(newposition-1)], collapse="|"),
                        as.character(results$Adducts[r]),
                        paste(adduc[newposition:length(adduc)],
                              collapse = "|"),
                        sep="|")
              }
            }
          }
        }
      }
    }
  }
  return(rawPeaks)
}


