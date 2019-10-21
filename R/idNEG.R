# idNEG
#' Lipids annotation for ESI-
#'
#' Lipids annotation based on fragmentation patterns for LC-MS/MS all-ions data
#' acquired in negative mode. This function compiles all functions writen for
#' ESI- annotations.
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
#' @param coelCutoff coelution score threshold between parent and fragment ions.
#' Only applied if rawData info is supplied. By default, 0.8.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be supplied. See \link{createLipidDB}
#' and \link{assignDB}.
#'
#' @return The output is a list with 2 elements: 1) a data frame that shows: ID,
#' class of lipid, CDB (total number of carbons and double bounds), FA
#' composition (specific chains composition if it has been confirmed), mz, RT
#' (in seconds), I (intensity, which comes directly from de input), Adducts,
#' ppm (m.z error), confidenceLevel (Subclass, FA level, where chains are known
#' but not their positions, or FA position level) and PFCS (parent-fragment
#' coelution score mean of all fragments used for the identification); and
#' 2) the original MS1 peaklist with the annotations on it.
#'
#' @examples
#' \donttest{
#' library(LipidMSdata)
#' idNEG(MS1 = MS1_neg, MSMS1 = MSMS1_neg, MSMS2 = MSMS2_neg)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idNEG <- function(MS1, MSMS1, MSMS2, ppm_precursor = 10,
                  ppm_products = 10, rttol = 10, coelCutoff = 0.8, dbs){

  if (missing(dbs)){
    dbs <- assignDB()
  }
  if (missing(MSMS2)){MSMS2 <- NULL}
  results <- vector()
  print("Searching for FA...")
  results <- rbind(results, idFAneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for FAHFA...")
  results <- rbind(results, idFAHFAneg(MS1 = MS1,
                                       MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                       ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                       dbs = dbs)$results)
  print("Searching for LPC...")
  results <- rbind(results, idLPCneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for LPE...")
  results <- rbind(results, idLPEneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for LPG...")
  results <- rbind(results, idLPGneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for LPI...")
  results <- rbind(results, idLPIneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for LPS...")
  results <- rbind(results, idLPSneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for PC...")
  results <- rbind(results, idPCneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for PE...")
  results <- rbind(results, idPEneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for PG...")
  results <- rbind(results, idPGneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for PI...")
  results <- rbind(results, idPIneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for PS...")
  results <- rbind(results, idPSneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for Sph...")
  results <- rbind(results, idSphneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for SphP...")
  results <- rbind(results, idSphPneg(MS1 = MS1,
                                      MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                      ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                      dbs = dbs)$results)
  print("Searching for Cer...")
  results <- rbind(results, idCerneg(MS1 = MS1,
                                     MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                     ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                     dbs = dbs)$results)
  print("Searching for CL...")
  results <- rbind(results, idCLneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Searching for Bile Acids...")
  results <- rbind(results, idBAneg(MS1 = MS1,
                                    MSMS1 = MSMS1, MSMS2 = MSMS2, ppm_precursor = ppm_precursor,
                                    ppm_products = ppm_products, rttol = rttol, coelCutoff = coelCutoff,
                                    dbs = dbs)$results)
  print("Preparing output...")
  if (nrow(results) > 0){
    annotatedPeaklist <- crossTables(MS1, results, ppm_precursor, rttol, dbs)
  } else {
    annotatedPeaklist <- data.frame()
  }
  return(list(results = results, annotatedPeaklist = annotatedPeaklist))
}
