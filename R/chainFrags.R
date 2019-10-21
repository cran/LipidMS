# chainFrags
#' Search of chain specific fragments
#'
#' Search of specific fragments that inform about the chains structure.
#'
#' @param coelfrags coeluting fragments for each candidate. Output of
#' \link{coelutingFrags}.
#' @param chainfrags character vector containing the fragmentation rules for
#' the chain fragments. If it is an empty vector, chains will be calculated based
#' on the difference between the precursor and the other chain. See details.
#' @param ppm m/z tolerance in ppm.
#' @param candidates candidates data frame. If any chain needs to be calculated
#' based on the difference between the precursor and the other chain, this
#' argument will be required. Output of \link{chainFrags}.
#' @param f known chains. If any chain needs to be calculated
#' based on the difference between the precursor and the other chain, this
#' argument will be required. Output of \link{chainFrags}.
#' @param dbs list of data bases required for the annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be changed. If data bases have
#' been customized using \link{createLipidDB}, they also have to be modified
#' here.
#'
#' @details The chainfrags argument must contain the fragmentation rules
#' which inform about the chains structure. For example, in the case of PG
#' subclass, the chain in sn1 position is identified by the lysoPG as M-H
#' resulting from the loss of the FA chain of sn2; and the chain in sn2 position
#' is identified as the free FA chain as M-H. These two fragments need to be
#' searched in two different steps: in the fist step we will look for lysoPGs
#' coeluting with the precursor using chainfrags = c("lysopg_M-H");
#' then, we will look for FA chains using chainfrags = c("fa_M-H"). This
#' information can be combined later using \link{combineChains} function.
#'
#' To indicate the fragments to be searched, the class of lipid is writen
#' using the same names as the LipidMS databases without the "db" at the end
#' (i.e. pa, dg, lysopa, mg, CE, etc.), and the adduct has to be indicated as
#' it appears in the adductsTable, both parts separated by "_". In case some
#' chain needs to be searched based on a neutral loss, this can be defined using
#' "NL-" prefix, followed by the database and adduct. If this neutral loss
#' is employed to find the remaining chain, "cbdiff-" prefix allows to calculate
#' the difference in carbons and doubles bounds between the precursor and the
#' building block found. For example, "cbdiff-dg_M+H-H2O" will look for DG as
#' M+H-H2O and then, it will return the difference between their number of
#' carbons and double bounds and the ones from the precursor. Otherwise,
#' "NL-mg_M+H-H2O" will look for fragments coming from the loss of MGs.
#'
#' In case these fragments identified as losses from the precursors are going to
#' be employed for the intensity rules, this same prefix has to be added.
#'
#' If a chain is calculated based on the difference of total number of
#' carbons and double bounds between the precursor and a previously searched chain,
#' \code{chainfrags} argument must be a character vector c("") and candidates
#' data frame and chain fragments list must be provided.
#'
#' @return List of data frames with the chain fragments found.
#'
#' @examples
#' \donttest{
#' library(LipidMSdata)
#' dbs <- assignDB()
#'
#' candidates <- findCandidates(MS1 = MS1_neg$peaklist,
#' db = dbs$pgdb, ppm = 10, rt = c(0, 2000), adducts = c("M-H"),
#' rttol = 10, rawData = MS1_neg$rawScans, coelCutoff = 0.8)
#'
#' MSMS <- rbind(MSMS1_neg$peaklist, MSMS2_neg$peaklist)
#' rawData <- rbind(MS1_neg$rawScans, MSMS1_neg$rawScans, MSMS2_neg$rawScans)
#' coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol = 10, rawData = rawData,
#' coelCutoff = 0.8)
#'
#' sn1 <- chainFrags(coelfrags, chainfrags = c("lysopg_M-H"), ppm = 10,
#' dbs = dbs)
#' sn2 <- chainFrags(coelfrags, chainfrags = c("fa_M-H"), ppm = 10, dbs = dbs)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
chainFrags = function (coelfrags, chainfrags, ppm = 10, candidates,
                       f = NULL, dbs){
  adductsTable <- dbs$adductsTable
  matches <- list()
  for (c in 1:nrow(candidates)) {
    fragments <- vector()
    for (i in 1:length(chainfrags)) {
      if (chainfrags[i] == ""){
        if (nrow(f[[c]]) > 0) {
          diffs <- diffcb(candidates$cb[c], f[[c]]$cb,
                          dbs$fadb)
          fmatch <- unique(data.frame(cb = diffs, m.z = 0, RT = f[[c]]$RT,
                               int = 0, peakID = "",
                               coelScore = NA,
                               db = "cbdiff", adduct = "",
                               stringsAsFactors = F))
          fmatch <- fmatch[fmatch$cb != "", ]
          fragments <- rbind(fragments, fmatch)
        }
      }
      if (grepl("cbdiff-", chainfrags[i])) {
        rule <- unlist(strsplit(chainfrags[i], "cbdiff-"))[2]
        db <- dbs[[paste(unlist(strsplit(rule, "_"))[1],
                         "db", sep = "")]]
        ad <- adductsTable[adductsTable$adduct ==
                                unlist(strsplit(rule, "_"))[2], ]
        mdiff <- ad$mdiff
        charge <- ad$charge
        n <- ad$n
        fmatch <- frags(coelfrags[[c]], ppm=ppm, db, mdiff, charge, n)
        if (nrow(fmatch) > 0) {
          fragments <-
            rbind(fragments,
                  data.frame(fmatch, db = paste(unlist(strsplit(chainfrags[i],
                            "_"))[1], sep = ""),
                            adduct = unlist(strsplit(chainfrags[i], "_"))[2]))
        }
      } else if (grepl("NL-", chainfrags[i])) {
        rule <- unlist(strsplit(chainfrags[i], "NL-"))[2]
        db <- dbs[[paste(unlist(strsplit(rule, "_"))[1],
                         "db", sep = "")]]
        db$Mass <- candidates$m.z[c] - db$Mass
        ad <- adductsTable[adductsTable$adduct ==
                                unlist(strsplit(rule, "_"))[2],]
        mdiff <- ad$mdiff
        charge <- ad$charge
        n <- ad$n
        fmatch <- frags(coelfrags[[c]], ppm, db, mdiff, charge, n)
        if (nrow(fmatch) > 0) {
          fragments <-
            rbind(fragments,
                  data.frame(fmatch,
                             db = paste(unlist(strsplit(chainfrags[i],
                                                        "_"))[1], sep = ""),
                             adduct = unlist(strsplit(chainfrags[i], "_"))[2]))
        }
      } else {
        rule <- chainfrags[i]
        db <- dbs[[paste(unlist(strsplit(rule, "_"))[1],
                         "db", sep = "")]]
        if (grepl("Mn", unlist(strsplit(chainfrags[i],
                                        "_"))[2])) {
          mdiff <- as.numeric(unlist(strsplit(unlist(strsplit(chainfrags[i],
                                                              "_"))[2], "Mn")))[2]
          charge <- 1
          n <- 1
        } else {
          ad <- adductsTable[adductsTable$adduct ==
                               unlist(strsplit(chainfrags[i], "_"))[2],]
          mdiff <- ad$mdiff
          charge <- ad$charge
          n <- ad$n
        }
        fmatch <- frags(coelfrags[[c]], ppm, db, mdiff, charge, n)
        if (nrow(fmatch) > 0) {
          fragments <- rbind(fragments,
                             data.frame(fmatch,
                                  db = unlist(strsplit(chainfrags[i], "_"))[1],
                                  adduct = unlist(strsplit(chainfrags[i],
                                                                 "_"))[2]))
        }
      }
    }
    if (class(fragments) == "data.frame") {
      matches[[c]] = unique(fragments)
    } else {
      matches[[c]] = data.frame()
    }
  }
  if (any(grepl("cbdiff", chainfrags))) {
    new <- list()
    for (c in 1:nrow(candidates)) {
      if (nrow(matches[[c]]) > 0) {
        fragments <- matches[[c]][grepl("cbdiff", matches[[c]]$db),
                                  ]
        other <- matches[[c]][!grepl("cbdiff", matches[[c]]$db),]
        db <- diffcb(candidates$cb[c], fragments$cb, dbs$fadb)
        fragments$cb <- db
        new[[c]] <- rbind(other, fragments[fragments$cb != "", ])
      } else {
        new[[c]] <- data.frame()
      }
    }
    matches <- new
  }
  return(matches)
}







