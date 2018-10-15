# idCLneg
#' Cardiolipines (CL) annotation for ESI-
#'
#' CL identification based on fragmentation patterns for LC-MS/MS
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
#' @param adducts expected adducts for CL in ESI-. Adducts allowed can
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
#' the chain fragments in sn1 position. See \link{chainFrags} for details.
#' @param chainfrags_sn2 character vector containing the fragmentation rules for
#' the chain fragments in sn2 position. See \link{chainFrags} for details.
#' @param chainfrags_sn3 character vector containing the fragmentation rules for
#' the chain fragments in sn3 position. See \link{chainFrags} for details.
#' @param chainfrags_sn4 character vector containing the fragmentation rules for
#' the chain fragments in sn4 position. See \link{chainFrags} for details.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}. If some intensity rules should be employed to
#' identify the chains position but they are't known yet, use "Unknown". If it
#' isn't required, leave an empty vector.
#' @param rates character vector with the expected rates between fragments given
#' as a string (i.e. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return List with CL annotations (results) and some additional information
#' (class fragments and chain fragments).
#'
#' @details \code{idCLneg} function involves 5 steps. 1) FullMS-based
#' identification of candidate CL as M-H or M-2H. 2) Search of CL class fragments:
#' no class fragments are searched by defaults as they use to have bad coelution
#' scores. 3) Search of specific fragments that inform about chain composition
#' at sn1 (lysoPA as M-H-H2O), sn2 (lysoPA as M-H-H2O), sn3 (lysoPA as M-H-H2O)
#' and sn4 (lysoPA as M-H-H2O). 4) Look for possible chains structure based on
#' the combination of chain fragments. 5) Check intensity rules to confirm
#' chains position. For CL there are no intensity rules by default.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level) and PFCS (parent-fragment coelution score mean of all fragments used
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
#' idCL(MS1 = LipidMS::MS1_neg, MSMS1 = LipidMS:MSMS1_neg,
#' MSMS2 = LipidMS::MSMS2_neg, coelCutoff = 0)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idCLneg <- function(MS1, MSMS1, MSMS2, ppm_precursor = 5,
                    ppm_products = 10, rttol = 5, rt,
                    adducts = c("M-H", "M+Na-2H"),
                    clfrags = c(),
                    clrequired = c(),
                    ftype = c(),
                    chainfrags_sn1 = c("lysopa_M-H-H2O"),
                    chainfrags_sn2 = c("lysopa_M-H-H2O"),
                    chainfrags_sn3 = c("lysopa_M-H-H2O"),
                    chainfrags_sn4 = c("lysopa_M-H-H2O"),
                    intrules = c("Unknown"),
                    rates = c(),
                    intrequired = c(),
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
  candidates <- findCandidates(MS1, dbs$cldb, ppm = ppm_precursor, rt = rt,
                               adducts = adducts, rttol = rttol, dbs = dbs,
                               rawData = rawData, coelCutoff = coelCutoff)

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
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    sn3 <- chainFrags(coelfrags, chainfrags_sn3, ppm_products, candidates, sn1,
                      dbs)
    sn4 <- chainFrags(coelfrags, chainfrags_sn4, ppm_products, candidates, sn1,
                      dbs)

    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=4, sn1, sn2, sn3, sn4)

    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=4,
                                   chainsComb)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 4, class="CL")


    if (length(clfrags) > 0){
      classfragments <- classConf$presence
      colnames(classfragments) <- clfrags
    } else {
      classfragments <- data.frame()
    }
    if (length(chainfrags_sn1) > 0){
      chainfragments <- mapply(rbind, sn1, sn2, sn3, sn4, SIMPLIFY=FALSE)
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
