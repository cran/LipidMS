# dataProcessing
#' Process mzXML files: peakpicking and deisotoping
#'
#' Process mzXML files: peak-picking using enviPick and deisotoping
#' using CAMERA.
#'
#' @param file path of the mzXML input file
#' @param msLevel numeric value indicating if data belongs to level 1 (fullMS)
#' or level 2 (MS/MS)
#' @param polarity character value: negative or positive
#' @param dmzgap enviPick parameter. 50 by default.
#' @param drtgap enviPick parameter. 25 by default.
#' @param ppm enviPick parameter. TRUE by default.
#' @param minpeak enviPick parameter. Optional. By default, 5 when msLevel = 1
#' and 4 when msLevel = 2.
#' @param maxint enviPick parameter. 1E9 by default.
#' @param dmzdens enviPick parameter. Optional. By default, 15 when msLevel = 1
#' and 20 when msLevel = 2.
#' @param drtdens enviPick parameter. Optional. 20 by default.
#' @param merged enviPick parameter. FALSE by default.
#' @param drtsmall enviPick parameter. Optional. By default, 100 when
#' msLevel = 1 and 20 when msLevel = 2.
#' @param drtfill enviPick parameter. 20 by default.
#' @param drttotal enviPick parameter. 200 by default.
#' @param recurs enviPick parameter. 4 by default.
#' @param weight enviPick parameter. Optional. By default, 1 when msLevel = 1
#' and 2 when msLevel = 2.
#' @param SB enviPick parameter. Optional. By default, 3 when msLevel = 1
#' and 2 when msLevel = 2.
#' @param SN enviPick parameter. 2 by default.
#' @param minint enviPick parameter. Optional. By default, 1000 when msLevel = 1
#' and 100 when msLevel = 2.
#' @param ended enviPick parameter. 2 by default.
#' @param progbar enviPick parameter. FALSE by default.
#' @param from enviPick parameter. FALSE by default.
#' @param to enviPick parameter. FALSE by default.
#'
#' @return Data frame with 3 columns (m.z, RT and int) containing the peaklist.
#'
#' @details This function executes 2 steps: 1) peak-picking using enviPick
#' package and 2) it removes isotopes using CAMERA. If msLevel = 1, only ions
#' with more than 1 isotope are kept.
#'
#' @examples
#'
#' \donttest{
#' dataProcessing("input_file.mzXML", msLevel = 1, polarity = "positive")
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
dataProcessing <- function(file, msLevel, polarity, dmzgap = 50, drtgap = 25,
                           ppm = TRUE, minpeak, maxint = 1E9, dmzdens, drtdens = 20,
                           merged = FALSE, drtsmall, drtfill = 20, drttotal = 200,
                           recurs = 4, weight, SB, SN = 2, minint, ended = 2,
                           progbar = FALSE, from = FALSE, to = FALSE){
  if (!msLevel %in% c(1, 2)){
    stop("msLevel must be 1 or 2")
  }
  if (!polarity %in% c("positive", "negative")){
    stop("polarity must be positive or negative")
  }
  if (file.exists(file)){
    if(missing(minpeak)){
      if (msLevel == 1){
        minpeak <- 5
      } else {
        minpeak <- 4
      }
    }
    if(missing(dmzdens)){
      if (msLevel == 1){
        dmzdens <- 15
      } else {
        dmzdens <- 20
      }
    }
    if(missing(minint)){
      if (msLevel == 1){
        minint <- 1000
      } else {
        minint <- 100
      }
    }
    if(missing(weight)){
      if (msLevel == 1){
        weight <- 2
      } else {
        weight <- 1
      }
    }
    if(missing(SB)){
      if (msLevel == 1){
        SB <- 3
      } else {
        SB <- 2
      }
    }
    if(missing(drtsmall)){
      if (msLevel == 1){
        drtsmall <- 100
      } else {
        drtsmall <- 20
      }
    }
    # peakPicking
    MSlist <- enviPick::readMSdata(file, MSlevel=1)
    MSlist <- enviPick::mzagglom(MSlist, dmzgap = dmzgap, ppm = ppm,
      drtgap = drtgap, minpeak = minpeak, maxint = maxint)
    MSlist <- enviPick::mzclust(MSlist, dmzdens = dmzdens, ppm = ppm,
      drtdens = drtdens, minpeak = minpeak, merged = merged)
    MSlist <- enviPick::mzpick(MSlist, minpeak = minpeak, drtsmall = drtsmall,
      drtfill = drtfill, drttotal = drttotal, recurs = recurs, weight = weight,
      SB = SB, SN = SN, minint = minint, maxint = maxint, ended = ended,
      progbar = progbar, from = from, to = to)
  # deisotope
  ## search isotopes
  sample <- makexcmsSet(MSlist, files=file, pD="data",
    polarity=polarity)
  sa <- CAMERA::xsAnnotate(sample, polarity=polarity)
  saF <- CAMERA::groupFWHM(sa, perfwhm = 0.8)
  saFI <- CAMERA::findIsotopes(saF, mzabs = 0.01, ppm = dmzdens)
  peaklist <- CAMERA::getPeaklist(saFI)
  ## remove isotopes
  M <- grepl("\\[M\\]\\+|\\[M\\]-", peaklist$isotopes)
  if (msLevel != 1){
    M <- grepl("\\[M\\]\\+|\\[M\\]-|^$", peaklist$isotopes)
  }
  peaks <- peaklist[M, c("mz", "rt", "into")]
  colnames(peaks) <- c("m.z", "RT", "int")
  return(peaks)
  } else {
    stop("The file doesn't exist!")
  }
}
