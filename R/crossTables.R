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
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return Data frame with 6 columns: m.z, RT, int, LipidMS_id, adduct and
#' confidence level for the annotation. When multiple IDs are proposed for the
#' same feature, they are sorted based on the annotation level.
#'
#' @examples
#' \donttest{
#' results <- idNEG(MS1 = MS1_neg, MSMS1 = MSMS1_neg, MSMS2 = MSMS2_neg)
#' crossTables(MS1_neg$peaklist, results = results$results,
#' ppm = 10, rttol = 10)}
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
crossTables <- function(MS1, results, ppm=10, rttol=10, dbs){

  # load dbs
  if (missing(dbs)){
    dbs <- assignDB()
  }
  #reorder MS data
  if (class(MS1) == "list" & length(MS1) == 2){
    if (!all(c("peaklist", "rawScans") %in% names(MS1))){
      stop("MS1, MSMS1 and MSMS2 (if supplied) lists should have two elements
           named as peaklist and rawScans")
    }
    if(!all(c("m.z", "RT", "int", "peakID") %in% colnames(MS1$peaklist))){
      stop("peaklist element of MS1, MSMS1 and MSMS2 needs to have at least 4
           columns: m.z, RT, int and peakID")
    }
    if(!all(c("m.z", "RT", "int", "peakID", "Scan") %in% colnames(MS1$rawScans))){
      stop("rawScans element of MS1, MSMS1 and MSMS2 needs to have at least 5
           columns: m.z, RT, int, peakID and Scan")
    }
    MS1 <- MS1$peaklist
    MS1$peakID <- as.vector(paste(MS1$peakID, "MS1", sep = "_"))
  }
  results <- results[order(results$CDB,
                           factor(results$confidenceLevel,
                                  levels= c("FA position", "FA",
                                            "Subclass", "MSMS",
                                            "MS-only")), -results$PFCS),]
  adductsTable <- dbs$adductsTable
  confLevels <- LipidMS::confLevels
  rawPeaks <- data.frame(MS1, LipidMS_id = "", Adduct = "",
                         confidenceLevel = "",
                         PFCS = "", stringsAsFactors = F)
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
        if (results$peakID[r] != "" && m == 1){
          rows <- which(MS1$peakID == results$peakID[r])
          matches <- c(1:length(rows))
        } else {
          if (sum(abs(rawPeaks$RT- as.numeric(results$RT[r])) < rttol) > 0){
            rows <- which(abs(rawPeaks$RT - as.numeric(results$RT[r])) < rttol)
            matches <- as.numeric(unlist(sapply(mzs[m], mzMatch,
                                                rawPeaks$m.z[rows], ppm)))
            matches <- matches[seq(1, length(matches), 2)]
          } else {
            matches <- vector()
          }
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
              rawPeaks$PFCS[rows[matches[i]]] <-
                as.character(results$PFCS[r])
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
                rawPeaks$PFCS[rows[matches[i]]] <-
                  paste(rawPeaks$PFCS[rows[matches[i]]],
                        as.character(results$PFCS[r]), sep="|")
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
                rawPeaks$PFCS[rows[matches[i]]] <-
                  paste(as.character(results$PFCS[r]),
                        rawPeaks$PFCS[rows[matches[i]]], sep="|")
              } else {
                ids <- unlist(strsplit(as.character(rawPeaks$LipidMS_id[rows[matches[i]]]),
                                       "\\|"))
                conflev <- unlist(strsplit(as.character(rawPeaks$confidenceLevel[rows[matches[i]]]),
                                           "\\|"))
                adduc <- unlist(strsplit(as.character(rawPeaks$Adduct[rows[matches[i]]]),
                                         "\\|"))
                score <- unlist(strsplit(as.character(rawPeaks$PFCS[rows[matches[i]]]),
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
                rawPeaks$PFCS[rows[matches[i]]] <-
                  paste(paste(score[1:(newposition-1)], collapse="|"),
                        as.character(results$PFCS[r]),
                        paste(score[newposition:length(score)],
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


