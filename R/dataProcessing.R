# dataProcessing
#' Process mzXML files: peakpicking and deisotoping
#'
#' Process mzXML files: peak-picking using enviPick and deisotoping
#' using an adaptation of the CAMERA algorithm.
#'
#' @param file path of the mzXML input file.
#' @param mslevel numeric value indicating if data belongs to level 1 (fullMS)
#' or level 2 (MS/MS).
#' @param polarity character value: negative or positive.
#' @param dmzgap enviPick parameter. 50 by default.
#' @param drtgap enviPick parameter. 25 by default.
#' @param ppm logical value. TRUE if dmzdens was set in ppm and FALSE if it was
#' in as an absolute value. TRUE by default.
#' @param minpeak minimum number of measurements required within the RT window
#' of drtsmall. Optional. By default, 5 when mslevel = 1 and 4 when mslevel = 2.
#' @param maxint EIC cluster with measurements above this intensity are kept,
#' even if they do not fulfill minpeak. 1E9 by default.
#' @param dmzdens maximum measurement deviation (+/-) of m/z from its mean
#' within each EIC. Optional. By default, 15 when mslevel = 1 and 30 when
#' mslevel = 2.
#' @param drtdens RT tolerance for clustering. Optional. 20 by default.
#' @param merged merge EIC cluster of comparable m/z. Logical. FALSE by default.
#' @param drtsmall peak definition - RT window of a peak. Optional. By default,
#' 100 when mslevel = 1 and 30 when mslevel = 2.
#' @param drtfill maximum RT gap length to be filled. 5 by default.
#' @param drttotal maximum RT length of a single peak. 100 by default.
#' @param recurs maximum number of peaks within one EIC. 3 by default.
#' @param weight weight for assigning measurements to a peak. Optional.
#' By default, 1 when mslevel = 1 and 2 when mslevel = 2.
#' @param SB signal-to-base ratio. Optional. By default, 3 when mslevel = 1
#' and 2 when mslevel = 2.
#' @param SN signal-to-noise ratio. 2 by default.
#' @param minint minimum intensity of a peakr. Optional. By default, 1000 when
#' mslevel = 1 and 100 when mslevel = 2.
#' @param ended within the peak detection recursion set by argument recurs,
#' how often can a peak detection fail to end the recursion?. 2 by default.
#' @param removeIsotopes logical. If TRUE, only isotopes identified as M+0, are
#' kept when mslevel = 1, and M+0 or unknown when mslevel = 2. TRUE by default.
#' If FALSE, an additional column is added to the peak list to inform about
#' isotopes.
#' @param rttolIso numeric. Time windows for isotope matching.
#' @param ppmIso numeric. Mass tolerance for isotope matching.
#'
#' @return List with two data frames: peaklist, with 4 columns (m.z, RT, int,
#' and peakID) and rawScan, with all the scans information in 5 columns (m.z,
#' RT, int, peakID and Scan). PeakID columns links both data frames: extracted
#' peaks and raw data. The Scan column indicates the scan number (order) to which
#' each row of the rawScans data frame belong.
#'
#' @details This function executes 2 steps: 1) peak-picking using enviPick
#' package and 2) it searches isotopes using an adaptation of the CAMERA
#' algorithm. If mslevel = 1 and remove isotopes is set as TRUE, only ions with
#' more than 1 isotope are kept.
#'
#' @examples
#' \donttest{
#' dataProcessing("input_file.mzXML", mslevel = 1, polarity = "positive")
#' }
#'
#' @references https://cran.r-project.org/web/packages/enviPick/index.html
#'
#' Kuhl C, Tautenhahn R, Boettcher C, Larson TR ans Neumann S (2012). "CAMERA:
#' an integrated strategy for compound spectra extraction and annotation of
#' liquid chromatography-mass spectrometry data sets." Analytical Chemistry, 84,
#' pp. 283-289. htto://pubs.acs.org/doi/abs/10.1021/ac202450g.
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
dataProcessing <- function(file, mslevel, polarity, dmzgap = 50, drtgap = 25,
                           ppm = TRUE, minpeak, maxint = 1E9, dmzdens, drtdens = 20,
                           merged = FALSE, drtsmall, drtfill = 5, drttotal = 100,
                           recurs = 4, weight, SB, SN = 2, minint, ended = 2,
                           removeIsotopes = TRUE, rttolIso = 2, ppmIso = 20){
  if (!mslevel %in% c(1, 2)){
    stop("mslevel must be 1 or 2")
  }
  if (!polarity %in% c("positive", "negative")){
    stop("polarity must be positive or negative")
  }
  if (file.exists(file)){
    if(missing(minpeak)){
      if (mslevel == 1){
        minpeak <- 5
      } else {
        minpeak <- 3
      }
    }
    if(missing(dmzdens)){
      if (mslevel == 1){
        dmzdens <- 15
      } else {
        dmzdens <- 30
      }
    }
    if(missing(minint)){
      if (mslevel == 1){
        minint <- 1000
      } else {
        minint <- 100
      }
    }
    if(missing(weight)){
      if (mslevel == 1){
        weight <- 2
      } else {
        weight <- 1
      }
    }
    if(missing(SB)){
      if (mslevel == 1){
        SB <- 3
      } else {
        SB <- 2
      }
    }
    if(missing(drtsmall)){
      if (mslevel == 1){
        drtsmall <- 100
      } else {
        drtsmall <- 30
      }
    }
    # peakPicking
    cat("\n Searching for features...")
    MSlist <- enviPick::readMSdata(file, MSlevel=1)
    MSlist <- enviPick::mzagglom(MSlist, dmzgap = dmzgap, ppm = ppm,
                                 drtgap = drtgap, minpeak = minpeak,
                                 maxint = maxint)
    MSlist <- enviPick::mzclust(MSlist, dmzdens = dmzdens, ppm = ppm,
                                drtdens = drtdens, minpeak = minpeak,
                                merged = merged)
    MSlist <- enviPick::mzpick(MSlist, minpeak = minpeak, drtsmall = drtsmall,
                               drtfill = drtfill, drttotal = drttotal,
                               recurs = recurs, weight = weight,
                               SB = SB, SN = SN, minint = minint, maxint = maxint,
                               ended = ended)
    cat("OK")
    peaklist <- data.frame(MSlist$Peaklist[,c("m/z", "RT", "sum_int", "peak_ID")],
                           stringsAsFactors = F)
    colnames(peaklist) <- c("m.z", "RT", "int", "peakID")
    rawScans <- data.frame(MSlist$Scans[[2]][,c("m/z", "RT", "intensity",
                                                "peakID")], stringsAsFactors = F)
    rawScans$Scan <- as.numeric(as.factor(rawScans$RT))
    colnames(rawScans) <- c("m.z", "RT", "int", "peakID", "Scan")

    cat("\n Annotating isotopes...")
    peaklistIso <- annotateIsotopes(peaklist, rawScans, ppm = ppmIso,
                                    rttol = rttolIso)
    cat("OK")

    if (removeIsotopes){
      cat("\n Removing isotopes...")
      if (mslevel == 1){
        peaklistIso <- peaklistIso[peaklistIso[,"isotope"] %in% c("[M+0]"),]
      } else {
        peaklistIso <- peaklistIso[peaklistIso[,"isotope"] %in% c("", "[M+0]"),]
      }

      cat("OK")
    }
    cat("\n")

    rawScans <- rawScans[rawScans[,"peakID"] %in% peaklistIso[,"peakID"],]

    return(list(peaklist = peaklistIso[,c("m.z", "RT", "int", "peakID")],
                rawScans = rawScans))
  } else {
    stop("The file doesn't exist!")
  }
}
