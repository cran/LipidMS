# idTGpos
#' Triacylglycerols (TG) annotation for ESI+
#'
#' TG identification based on fragmentation patterns for LC-MS/MS
#' AIF data acquired in positive mode.
#'
#' @param MS1 data frame cointaining all peaks from the full MS function. It
#' must have three columns: m.z, RT (in seconds) and int (intensity).
#' @param MSMS1 data frame cointaining all peaks from the low energy function.
#' It must have three columns: m.z, RT and int.
#' @param MSMS2 data frame cointaining all peaks from the high energy
#' function if it is the case. It must have three columns: m.z, RT and int.
#' Optional.
#' @param ppm_precursor mass tolerance for precursor ions. By default, 5 ppm.
#' @param ppm_products mass tolerance for product ions. By default, 10 ppm.
#' @param rttol total rt window for coelution between precursor and product
#' ions. By default, 3 seconds.
#' @param rt rt range where the function will look for candidates. By default,
#' it will search within all RT range in MS1.
#' @param adducts expected adducts for TG in ESI+. Adducts allowed can
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
#' the chain fragments in sn2 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn1 chains.
#' @param chainfrags_sn3 character vector containing the fragmentation rules for
#' the chain fragments in sn3 position. See \link{chainFrags} for details. If
#' empty, it will be estimated based on the difference between precursors and
#' sn2 chains.
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}. If some intensity rules should be employed to
#' identify the chains position but they are't known yet, use "Unknown". If it
#' isn't required, leave an empty vector.
#' @param rates character vector with the expected rates between fragments given
#' as a string (i.e. "3/1"). See \link{checkIntensityRules}.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be changed. If data bases have
#' been customized using \link{createLipidDB}, they also have to be modified
#' here.
#'
#' @return list with TG annotations (results) and some additional information
#' (class fragments and chain fragments).
#'
#' @details \code{idTGpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate TG as M+NH4 and M+Na (for long chain TGs
#' only those that are found as M+NH4 and M+Na are kept to avoid FP).
#' 2) Search of TG class fragments: there are no class fragment by default.
#' 3) Search of specific fragments that inform about the FA chains: DGs resulting
#' from the loss of FA chains as M+H-H2O. Only those with a relative intensity
#' to the precursor greater than 0.1 are kept to avoid FP. 4) Look for possible
#' chains structure based on the combination of chain fragments. 5) Check
#' intensity rules to confirm chains position. In the case of TG, FA position
#' level still needs some improvements.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (Subclass,
#' FA level, where chains are known but not their positions, or FA position
#' level).
#'
#' @note Isotopes should be removed before identification to avoid false
#' positives.
#' This function has been writen based on fragmentation patterns observed for
#' two different platforms (QTOF 6550 from Agilent and Sinapt G2-Si from Waters),
#' but it may need to be customized for other platforms or acquisition settings.
#'
#' @examples
#' \donttest{idTGpos(MS1 = LipidMS::mix_pos_fullMS, MSMS1 = LipidMS::mix_pos_Ce20,
#' MSMS2 = LipidMS::mix_pos_Ce40)}
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idTGpos <- function(MS1, MSMS1, MSMS2 = data.frame(), ppm_precursor = 5,
                    ppm_products = 10, rttol = 3,
                    rt = c(min(MS1$RT), max(MS1$RT)),
                    adducts = c("M+NH4", "M+Na"),
                    clfrags = c(),
                    clrequired = c(),
                    ftype = c(),
                    chainfrags_sn1 = c("cbdiff-dg_M+H-H2O"),
                    chainfrags_sn2 = c("cbdiff-dg_M+H-H2O"),
                    chainfrags_sn3 = c("cbdiff-dg_M+H-H2O"),
                    intrules = c("Unknown"),
                    rates = c(),
                    intrequired = c(T),
                    dbs = list(tgdb = LipidMS::tgdb,
                               fadb = LipidMS::fadb,
                               dgdb = LipidMS::dgdb,
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
  candidates <- findCandidates(MS1, dbs[["tgdb"]], ppm_precursor, rt, adducts,
                               rttol, dbs)
  if (all(c("M+NH4", "M+Na") %in% adducts)){
    candidates <- candidates[grepl("M\\+NH4;M\\+Na|M\\+Na;M\\+NH4",
                                candidates$adducts) | candidates$cb < "42:0",]
  }
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
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)
    sn3 <- chainFrags(coelfrags, chainfrags_sn3, ppm_products, candidates,
                      f = sn1, dbs)
    sn1 <- removeDGs(sn1, candidates)
    sn2 <- removeDGs(sn2, candidates)
    sn3 <- removeDGs(sn3, candidates)
    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=3, sn1, sn2, sn3)

    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=3,
                                   chainsComb, sn1, sn2, sn3)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 3, class="TG")

    if (length(clfrags) > 0){
      classfragments <- classConf$presence
      colnames(classfragments) <- clfrags
    } else {
      classfragments <- data.frame()
    }
    if (length(chainfrags_sn1) > 0){
      chainfragments <- mapply(rbind,sn1,sn2, sn3,SIMPLIFY=FALSE)
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
