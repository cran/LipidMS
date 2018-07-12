# idSphpos
#' Sphingoid bases (Sph) annotation for ESI-
#'
#' Sph identification based on fragmentation patterns for LC-MS/MS
#' AIF data acquired in positive mode.
#'
#' @param MS1 data frame cointaining all peaks from the full MS function. It
#' must have three columns: m.z, RT (in seconds) and int.
#' @param MSMS1 data frame cointaining all peaks from the first high energy
#' function. It must have three columns: m.z, RT (in seconds) and int.
#' @param MSMS2 data frame cointaining all peaks from the second high energy
#' function if it is the case. It must have three columns: m.z, RT (in seconds)
#' and int. Optional.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursors and product
#' ions. By default, 3 seconds.
#' @param rt rt window where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for Sph in ESI+. Adducts allowed can
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
#' @param dbs list of data bases required for the annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be changed. If data bases have
#' been customized using \link{createLipidDB}, they also have to be uploaded here.
#'
#' @return list with Sph annotations (results) and some additional information
#' (class fragments and chain fragments).
#'
#' @details \code{idSphpos} function involves 3 steps. 1) FullMS-based
#' identification of candidate Sph as M+H. 2) Search of Sph class fragments in
#' case some is defined. 3) Search of specific fragments that confirm chain
#' composition (neutral loss of 1 or 2 H2O molecules).
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (in this
#' case, as Sph only have one chain, only Subclass and FA level are possible).
#'
#' @note Isotopes should be removed before identification to avoid false
#' positives.
#' This function has been writen based on fragmentation patterns observed for
#' a Q-TOF 6550 from Agilent.
#'
#' @examples
#' idSphpos(MS1 = LipidMS::mix_pos_fullMS, MSMS1 = LipidMS::mix_pos_Ce20,
#' MSMS2 = LipidMS::mix_pos_Ce40)
#'
#' idSphpos(MS1 = LipidMS::serum_pos_fullMS, MSMS1 = LipidMS::serum_pos_Ce20,
#' MSMS2 = LipidMS::serum_pos_Ce40)
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idSphpos <- function(MS1, MSMS1, MSMS2 = data.frame(), ppm_precursor = 5,
                     ppm_products = 10, rttol = 3,
                     rt = c(min(MS1$RT), max(MS1$RT)),
                     adducts = c("M+H"),
                     clfrags = c(),
                     clrequired = c(),
                     ftype = c(),
                     chainfrags_sn1 = c("sph_M+H-H20", "sph_M+H-2H2O"),
                     dbs = list(sphdb = LipidMS::sphdb,
                                adductsTable = LipidMS::adductsTable)){

  if (!all(colnames(MS1) %in% c("m.z", "RT", "int")) |
      !all(colnames(MSMS1) %in% c("m.z", "RT", "int"))){
    stop("Peaklists (MS1, MSMS1 and MSMS2 if supplied, should have 3 columns
         with the following names: m.z, RT, int.")
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

  # candidates search
  candidates <- findCandidates(MS1, dbs[["sphdb"]], ppm_precursor, rt, adducts,
                               rttol, dbs)

  if (nrow(candidates) > 0){
    # isolation of coeluting fragments
    MSMS <- rbind(MSMS1, MSMS2)
    coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol)

    # check class fragments
    classConf <- checkClass(candidates, coelfrags, clfrags, ftype, clrequired,
                            ppm_products, dbs)

    # search chains fragments
    sn1 <- chainFrags(coelfrags, chainfrags_sn1, ppm_products, dbs = dbs,
                      candidates = candidates)

    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=1, sn1)

    # prepare output
    res <-organizeResults(candidates, clfrags, classConf, chainsComb,
                          intrules = c(), intConf = c(), nchains = 1,
                          class="Sph")

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
