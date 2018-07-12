# idPEpos
#' Phosphoethanolamines (PE) annotation for ESI+
#'
#' PE identification based on fragmentation patterns for LC-MS/MS
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
#' @param adducts expected adducts for PE in ESI+. Adducts allowed can
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
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
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
#' @return List with PE annotations (results) and some additional information
#' (class fragments and chain fragments).
#'
#' @details \code{idPEpos} function involves 5 steps. 1) FullMS-based
#' identification of candidate PE as M+H and M+Na. 2) Search of PE class
#' fragments: loss of head group (DG as M+H-H2O) coeluting with the precursor
#' ion. 3) Search of specific fragments that inform about chain composition at
#' sn1 (MG as M+H-H2O resulting from the loss of the FA chain at sn2 and
#' the head group or LPE as M+H-H2O resulting just from the loss of the FA chain)
#' and sn2 (FA or MG chain from sn2as M+H-H2O). 4) Look for possible chains
#' structure based on the combination of chain fragments. 5) Check intensity
#' rules to confirm chains position. LPE or MG from sn1 is at least 3 times more
#' intense than the ones from sn2.
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
#' idPEpos(MS1 = LipidMS::mix_pos_fullMS, MSMS1 = LipidMS::mix_pos_Ce20,
#' MSMS2 = LipidMS::mix_pos_Ce40)
#'
#' idPEpos(MS1 = LipidMS::serum_pos_fullMS, MSMS1 = LipidMS::serum_pos_Ce20,
#' MSMS2 = LipidMS::serum_pos_Ce40)
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idPEpos <- function(MS1, MSMS1, MSMS2 = data.frame(), ppm_precursor = 5,
                    ppm_products = 10, rttol = 3,
                    rt = c(min(MS1$RT), max(MS1$RT)),
                    adducts = c("M+Na", "M+H"),
                    clfrags = c("dg_M+H-H2O"),
                    clrequired = c(F),
                    ftype = c("BB"),
                    chainfrags_sn1 = c("lysope_M+H-H2O", "mg_M+H-H2O"),
                    chainfrags_sn2 = c("fa_M+H-H2O", "mg_M+H-H2O"),
                    intrules = c("mg_sn1/mg_sn1", "lysope_sn1/lysope_sn1"),
                    rates = c("3/1", "3/1"),
                    intrequired = c(T),
                    dbs = list(pedb = LipidMS::pedb,
                               dgdb = LipidMS::dgdb,
                               mgdb = LipidMS::mgdb,
                               lysopedb = LipidMS::lysopedb,
                               fadb = LipidMS::fadb,
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
  candidates <- findCandidates(MS1, dbs[["pedb"]], ppm_precursor, rt, adducts,
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
    sn2 <- chainFrags(coelfrags, chainfrags_sn2, ppm_products, candidates, sn1,
                      dbs)

    # combine chain fragments
    chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)

    # check chains position based on intensity ratios
    intConf <- checkIntensityRules(intrules, rates, intrequired, nchains=2,
                                   chainsComb, sn1, sn2)

    # prepare output
    res <- organizeResults(candidates, clfrags, classConf, chainsComb, intrules,
                           intConf, nchains = 2, class="PE")
    if (length(clfrags) > 0){
      classfragments <- classConf$presence
      colnames(classfragments) <- clfrags
    } else {
      classfragments <- data.frame()
    }
    if (length(chainfrags_sn1) > 0){
      chainfragments <- mapply(rbind,sn1,sn2,SIMPLIFY=FALSE)
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
