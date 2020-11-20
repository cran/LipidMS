# idLPSneg
#' Lysophosphoserines (LPS) annotation for ESI-
#'
#' LPS identification based on fragmentation patterns for LC-MS/MS
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
#' @param adducts expected adducts for LPS in ESI-. Adducts allowed can
#' be modified in adductsTable (dbs argument).
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See \link{checkClass} for details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See \link{checkClass} for details.
#' @param clrequired logical vector indicating if each class fragment is
#' required or not. If any of them is required, at least one of them must be
#' present within the coeluting fragments. See \link{checkClass} for details.
#' @param chainfrags_sn1 character vector containing the fragmentation rules for
#' the chain fragments. See \link{chainFrags} for details.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return List with LPS annotations (results) and some additional information
#' (class fragments and chain fragments).
#'
#' @details \code{idLPSneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate LPS as M-H and M+Na-2H. 2) Search of
#' LPS class fragments: neutral loss of 87.032 coeluting with the precursor ion.
#' 3) Search of specific fragments that confirm chain composition (FA as M-H).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as LPS only have one chain, only Subclass and FA level are possible)
#' and PFCS (parent-fragment coelution score mean of all fragments used
#' for the identification).
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
#' idLPSneg(MS1 = MS1_neg, MSMS1 = MSMS1_neg, MSMS2 = MSMS2_neg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idLPSneg <- function(MS1, MSMS1, MSMS2, ppm_precursor = 5,
                     ppm_products = 10, rttol = 3, rt,
                     adducts = c("M-H", "M+Na-2H"),
                     clfrags = c(87.032),
                     clrequired = c(F),
                     ftype = c("NL"),
                     chainfrags_sn1 = c("fa_M-H"),
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
  if (!any(c("m.z", "RT", "int") %in% colnames(MS1)) |
      !any(c("m.z", "RT", "int") %in% colnames(MSMS1))){
    stop("Peaklists (MS1, MSMS1 and MSMS2 if supplied, should have at least
         3 columns with the following names: m.z, RT, int.")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
  }
  if (length(clfrags) > 0){
    if (length(clfrags) != length(clrequired) | length(clfrags) !=
        length(ftype)){
      stop("clfrags, clrequired and ftype should have the same length")
    }
    if (!all(ftype %in% c("F", "NL", "BB"))){
      stop("ftype values allowed are: \"F\", \"NL\" or\"BB\"")
    }
    strfrag <- which(grepl("_", clfrags))
    if (length(strfrag) > 0){
      d <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 1))
      a <- unlist(lapply(strsplit(clfrags[strfrag], "_"), "[[", 2))
      if (!all(a %in% dbs[["adductsTable"]]$adduct)){
        stop("Adducts employed in clfrags also need to be at adductsTable.")
      }
      if (!all(paste(d, "db", sep="") %in% names(dbs))){
        stop("All required dbs must be supplied through dbs argument.")
      }
    }
  }
  if (missing(rt)){
    rt <- c(min(MS1$RT), max(MS1$RT))
  }

  rawData <- rbind(rawDataMS1, rawDataMSMS1, rawDataMSMS2)
  rawData <- rawData[!rawData$peakID %in% c("0_MS1", "0_MSMS1", "0_MSMS2"),]
  # candidates search
  candidates <- findCandidates(MS1, dbs$lysopsdb, ppm = ppm_precursor,
                                rt = rt, adducts = adducts, rttol = rttol,
                                dbs = dbs, rawData = rawData,
                               coelCutoff = coelCutoff)

  if (nrow(candidates) > 0){
    # isolation of coeluting fragments
    MSMS <- rbind(MSMS1, MSMS2)
    coelfrags <- coelutingFrags(candidates, MSMS, rttol, rawData,
                                coelCutoff = coelCutoff)

    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)

    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)

    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=1, sn1)

    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules = c(), rates = c(),
                                   intrequired = c(), nchains=1,
                                   chainsComb)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb,
                           intrules  = c(), intConf, nchains = 1, class="LPS")

    if (length(clfrags) > 0){
      classfragments <- classConf$presence
      colnames(classfragments) <- clfrags
    } else {
      classfragments <- data.frame()
    }
    if (length(chainfrags_sn1) > 0){
      chainfragments <- sn1
    } else {
      chainsfrags <- list()
    }

    return(list(results = res, candidates = candidates,
                classfragments = classfragments,
                chainfragments = chainfragments))
  } else {
    return(list(results = data.frame()))
  }
}
