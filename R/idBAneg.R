# idBAneg
#' Bile Acids (BA) annotation for ESI-
#'
#' BA identification based on fragmentation patterns for LC-MS/MS
#' AIF data acquired in negative mode.
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
#' @param adducts expected adducts for BA in ESI-. Adducts allowed can
#' be modified in the adducsTable (dbs argument).
#' @param conjfrag character vector containing the fragmentation rules for
#' the BA-conjugates. By default just taurine and glycine are considered,
#' but baconjdb can be modified to add more possible conjugates.
#' See \link{chainFrags} for details. It can also be an empty vector.
#' @param bafrag character vector containing the fragmentation rules for
#' other BA fragments. See \link{chainFrags} for details. It can be an empty
#' vector.
#' @param dbs list of data bases required for annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be changed. If data bases have
#' been customized using \link{createLipidDB}, they also have to be modified
#' here.
#'
#' @return List with BA annotations (results) and some additional information
#' (fragments).
#'
#' @details \code{idBAneg} function involves 3 steps. 1) FullMS-based
#' identification of candidate BA as M-H. 2) Search of BA-conjugate fragments if
#' required. 3) Search of fragments coming from the loss of H2O. If any
#' fragmentation rule has been defined, only those with MSMS support are kept.
#'
#' Results data frame shows: ID, class of lipid, CDB (total number
#' of carbons and double bounds), FA composition (specific chains composition if
#' it has been confirmed), mz, RT (in seconds), I (intensity, which comes
#' directly from de input), Adducts, ppm (m.z error), confidenceLevel (MS-only
#' if no rules are defined, or MSMS level if they are supported by fragments).
#'
#' @note Isotopes should be removed before identification to avoid false
#' positives.
#' This function has been writen based on fragmentation patterns observed for
#' two different platforms (QTOF 6550 from Agilent and Sinapt G2-Si from Waters),
#' but it may need to be customized for other platforms or acquisition settings.
#'
#' @examples
#' idBAneg(MS1 = LipidMS::mix_neg_fullMS, MSMS1 = LipidMS::mix_neg_Ce20,
#' MSMS2 = LipidMS::mix_neg_Ce40)
#'
#' idBAneg(MS1 = LipidMS::serum_neg_fullMS, MSMS1 = LipidMS::serum_neg_Ce20,
#' MSMS2 = LipidMS::serum_neg_Ce40)
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
idBAneg <- function(MS1, MSMS1, MSMS2 = data.frame(), ppm_precursor = 5,
                    ppm_products = 10, rttol = 3,
                    rt = c(min(MS1$RT), max(MS1$RT)),
                    adducts = c("M-H"),
                    conjfrag = c("baconj_M-H"),
                    bafrag = c("ba_M-H-H2O"),
                    dbs = list(badb = LipidMS::badb,
                               baconjdb = LipidMS::baconjdb,
                               adductsTable = LipidMS::adductsTable)){

  if (!all(colnames(MS1) %in% c("m.z", "RT", "int")) |
      !all(colnames(MSMS1) %in% c("m.z", "RT", "int"))){
    stop("Peaklists (MS1, MSMS1 and MSMS2 if supplied, should have 3 columns
         with the following names: m.z, RT, int.")
  }
  if (!all(adducts %in% dbs[["adductsTable"]]$adduct)){
    stop("Some adducts can't be found at the aductsTable. Add them.")
  }

  # candidates search
  candidates <- findCandidates(MS1, dbs[["badb"]], ppm_precursor, rt, adducts,
                               rttol, dbs)

  if (nrow(candidates) > 0){
    # isolation of coeluting fragments
    MSMS <- rbind(MSMS1, MSMS2)
    coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol)

    # search fragments
    check <- rep(list(vector()), nrow(candidates))
    if (length(conjfrag) > 0){
      conjugates <- chainFrags(coelfrags, conjfrag, ppm_products, dbs = dbs,
                               candidates = candidates)
      for (c in 1:nrow(candidates)){
        if (nrow(conjugates[[c]]) > 0){
          if (dbs[["badb"]][dbs[["badb"]]$total ==
                            candidates$cb[c], "conjugate"] %in%
              conjugates[[c]]$cb | dbs[["badb"]][dbs[["badb"]]$total ==
                                                 candidates$cb[c],
                                       "conjugate"] == ""){
            check[[c]] <- append(check[[c]], TRUE)
          } else {
            check[[c]] <- append(check[[c]], FALSE)
          }
        } else {
          check[[c]] <- append(check[[c]], FALSE)
        }
      }
    }
    if (length(bafrag) > 0){
      bas <- chainFrags(coelfrags, bafrag, ppm_products, dbs = dbs,
                        candidates = candidates)
      for (c in 1:nrow(candidates)){
        if(nrow(bas[[c]]) > 0){
          if (candidates$cb[c] %in% bas[[c]]$cb){
            check[[c]] <- append(check[[c]], TRUE)
          } else {
            check[[c]] <- append(check[[c]], FALSE)
          }
        } else {
          check[[c]] <- append(check[[c]], FALSE)
        }
      }
    }
    check <- matrix(unlist(check), nrow = nrow(candidates), byrow = T)
    colnames(check) <- c(conjfrag, bafrag)


    # keep only those supported by fragments
    if (length(bafrag) > 0 | length(conjfrag) > 0){
      candidates2keep <- candidates[unlist(apply(check, 1, sum)) > 0,]
    } else {
      candidates2keep <- candidates
    }

    # prepare output
    res <- organizeResults(candidates2keep, clfrags = c(), classConf = c(),
                           chainComb = c(), intrules = c(),
                           intConf = c(), nchains = 0, class="BA")
    if (length(bafrag) > 0 | length(conjfrag) > 0){
      res$confidenceLevel <- "MSMS"
    }
    return(list(results = res, candidates = candidates, fragments = check))
  } else {
    return(list(results = data.frame()))
  }
}
