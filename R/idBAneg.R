# idBAneg
#' Bile Acids (BA) annotation for ESI-
#'
#' BA identification based on fragmentation patterns for LC-MS/MS
#' AIF data acquired in negative mode.
#'
#' @param MS1 list with two data frames cointaining all peaks from the full MS
#' function ("peaklist" data frame) and the raw MS scans data ("rawScans" data
#' frame). They must have four columns: m.z, RT (in seconds), int (intensity)
#' and peakID (link between both data frames). "rawScans" data frame also needs
#' a extra column named "Scan", which indicates the scan order number. Output
#' of \link{dataProcessing} function. In case no coelution score needs to be
#' applied, this argument can be just the peaklist data frame.
#' @param MSMS1 list with two data frames cointaining all peaks from the high
#' energy function ("peaklist" data frame) and the raw MS scans data ("rawScans"
#' data frame). They must have four columns: m.z, RT (in seconds), int (intensity)
#' and peakID (link between both data frames). "rawScans" data frame also needs
#' a extra column named "Scan", which indicates the scan order number. Output
#' of \link{dataProcessing} function. In case no coelution score needs to be
#' applied, this argument can be just the peaklist data frame.
#' @param MSMS2 list with two data frames cointaining all peaks from a second high
#' energy function ("peaklist" data frame) and the raw MS scans data ("rawScans"
#' data frame). They must have four columns: m.z, RT (in seconds), int (intensity)
#' and peakID (link between both data frames). "rawScans" data frame also needs
#' a extra column named "Scan", which indicates the scan order number. Output
#' of \link{dataProcessing} function. In case no coelution score needs to be
#' applied, this argument can be just the peaklist data frame. Optional.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for BA in ESI-. Adducts allowed can
#' be modified in the adducsTable (dbs argument).
#' @param conjfrag character vector containing the fragmentation rules for
#' the BA-conjugates. By default just taurine and glycine are considered,
#' but baconjdb can be modified to add more possible conjugates.
#' See \link{chainFrags} for details. It can also be an empty vector.
#' @param bafrag character vector containing the fragmentation rules for
#' other BA fragments. See \link{chainFrags} for details. It can be an empty
#' vector.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return List with BA annotations (results) and some additional information
#' (fragments).
#'
#' @details \code{idBAneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate BA as M-H. 2) Search of BA-conjugate fragments if
#' required. 3) Search of fragments coming from the loss of H2O.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (MS-only
#' if no rules are defined, or Subclass level if they are supported by fragments)
#' and PFCS (parent-fragment coelution score mean of all fragments used for the
#' identification).
#'
#' @note Isotopes should be removed before identification to avoid false
#' positives.
#' This function has been writen based on fragmentation patterns observed for
#' two different platforms (QTOF 6550 from Agilent and Sinapt G2-Si from Waters),
#' but it may need to be customized for other platforms or acquisition settings.
#'
#' @examples
#' \donttest{
#' library(LipidMSdata)
#' idBAneg(MS1 = MS1_neg, MSMS1 = MSMS1_neg, MSMS2 = MSMS2_neg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idBAneg <- function(MS1, MSMS1, MSMS2, ppm_precursor = 5,
                    ppm_products = 10, rttol = 3, rt,
                    adducts = c("M-H"),
                    conjfrag = c("baconj_M-H"),
                    bafrag = c("ba_M-H-H2O", "ba_M-H-2H2O"),
                    coelCutoff = 0.8,
                    dbs){

  # load dbs
  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (missing(MSMS2) | is.null(MSMS2)){
    rawDataMSMS2 <- MSMS2 <- data.frame()
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
    rawDataMS1 <- MS1$rawScans
    rawDataMS1$peakID <- as.vector(paste(rawDataMS1$peakID, "MS1", sep = "_"))
    MS1 <- MS1$peaklist
    MS1$peakID <- as.vector(paste(MS1$peakID, "MS1", sep = "_"))
    rawDataMSMS1 <- MSMS1$rawScans
    rawDataMSMS1$peakID <- as.vector(paste(rawDataMSMS1$peakID, "MSMS1", sep = "_"))
    MSMS1 <- MSMS1$peaklist
    MSMS1$peakID <- as.vector(paste(MSMS1$peakID, "MSMS1", sep = "_"))
    if (!missing(MSMS2) & !is.data.frame(MSMS2)){
      rawDataMSMS2 <- MSMS2$rawScans
      rawDataMSMS2$peakID <- as.vector(paste(rawDataMSMS2$peakID, "MSMS2", sep = "_"))
      MSMS2 <- MSMS2$peaklist
      MSMS2$peakID <- as.vector(paste(MSMS2$peakID, "MSMS2", sep = "_"))
    }
  } else {
    rawDataMS1 <- rawDataMSMS1 <- rawDataMSMS2 <- data.frame()
    coelCutoff <- 0
  }
  if(!"peakID" %in% colnames(MS1)){
    MS1$peakID <- as.vector(rep("", nrow(MS1)))
    MSMS1$peakID <- as.vector(rep("", nrow(MSMS1)))
    if (nrow(MSMS2) != 0){
      MSMS2$peakID <- as.vector(rep("", nrow(MSMS2)))
    }
  }
  if (!(c("m.z", "RT", "int") %in% colnames(MS1)) ||
      !(c("m.z", "RT", "int") %in% colnames(MS1))){
    stop("Peaklists (MS1, MSMS1 and MSMS2 if supplied, should have at least
         3 columns with the following names: m.z, RT, int.")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
  }
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }

  rawData <- rbind(rawDataMS1, rawDataMSMS1, rawDataMSMS2)
  rawData <- rawData[!rawData$peakID %in% c("0_MS1", "0_MSMS1", "0_MSMS2"),]
  # candidates search
  candidates <- findCandidates(MS1, dbs$badb, ppm = ppm_precursor, rt = rt,
                                adducts = adducts, rttol = rttol, dbs = dbs,
                                rawData = rawData, coelCutoff = coelCutoff)


  if (nrow(candidates) > 0){
    # isolation of coeluting fragments
    MSMS <- rbind(MSMS1, MSMS2)
    coelfrags <- coelutingFrags(candidates, MSMS, rttol, rawData,
                                coelCutoff = coelCutoff)

    # search fragments
    check <- rep(list(vector()), nrow(candidates))
    conjugates <- rep(list(vector()), nrow(candidates))
    bas <- rep(list(vector()), nrow(candidates))
    if (length(conjfrag) > 0){
      conjugates <- chainFrags(coelfrags, conjfrag, ppm_products, dbs = dbs,
                               candidates = candidates)
      for (c in 1:nrow(candidates)){
        if (nrow(conjugates[[c]]) > 0){
          if (dbs$badb[dbs$badb$total ==
                       candidates$cb[c], "conjugate"] %in%
              conjugates[[c]]$cb | dbs$badb[dbs$badb$total ==
                                            candidates$cb[c],
                                            "conjugate"] == ""){
            check[[c]] <- append(check[[c]], TRUE)
          } else {
            check[[c]] <- append(check[[c]], FALSE)
          }
        } else if (dbs$badb[dbs$badb$total ==
                            candidates$cb[c],
                            "conjugate"] == ""){
          check[[c]] <- append(check[[c]], TRUE)
          conjugates[[c]] <- data.frame(0,0,0,0,0,0,0,0)
          colnames(conjugates[[c]]) <- c("cb", "m.z", "RT", "int",
                                         "peakID", "coelScore", "db", "adduct")
        } else {
          check[[c]] <- append(check[[c]], FALSE)
          conjugates[[c]] <- data.frame(0,0,0,0,0,0,0,0)
          colnames(conjugates[[c]]) <- c("cb", "m.z", "RT", "int",
                                         "peakID", "coelScore", "db", "adduct")
        }
      }
      classConf <- list()
    }
    if (length(bafrag) > 0){
      bas <- chainFrags(coelfrags, bafrag, ppm_products, dbs = dbs,
                        candidates = candidates)
      for (c in 1:nrow(candidates)){
        if(nrow(bas[[c]]) > 0){
          if (candidates$cb[c] %in% bas[[c]]$cb ||
              dbs$badb[dbs$badb$total == candidates$cb[c], "base"]
              %in% bas[[c]]$cb){
            check[[c]] <- append(check[[c]], TRUE)
          } else {
            check[[c]] <- append(check[[c]], FALSE)
          }
        } else {
          check[[c]] <- append(check[[c]], FALSE)
          bas[[c]] <- data.frame(0,0,0,0,0,NA,0,0)
          colnames(bas[[c]]) <- c("cb", "m.z", "RT", "int",
                                         "peakID", "coelScore", "db", "adduct")
        }
      }
    }
    check <- matrix(unlist(check), nrow = nrow(candidates), byrow = T)
    classConf$presence <- check
    classConf$passed <- unlist(apply(check, 1, sum)) > 0
    classConf$fragments <- Map(rbind, conjugates, bas)

    # prepare output
    if (length(bafrag) > 0 | length(conjfrag) > 0){
      clfrags <- c("fragment")
    } else {
      clfrags <- c()
    }
    res <- organizeResults(candidates, clfrags = clfrags, classConf = classConf,
                           chainsComb = c(), intrules = c(),
                           intConf = c(), nchains = 0, class="BA")
    return(list(results = res, candidates = candidates, classConf))
  } else {
    return(list(results = data.frame()))
  }
}
