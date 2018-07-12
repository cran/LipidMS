# organizeResults
#' Prepare output for LipidMS annotation functions
#'
#' Prepare a readable output for LipidMS identification functions.
#'
#' @param candidates candidates data frame. Output of \link{findCandidates}.
#' @param clfrags vector containing the expected fragments for a given lipid
#' class.
#' @param classConf output of \link{checkClass}
#' @param chainComb output of \link{combineChains}
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param intConf output of \link{checkIntensityRules}
#' @param nchains number of chains of the targeted lipid class.
#' @param class character value. Lipid class (i.e. PC, PE, DG, TG, etc.).
#'
#' @examples
#' \donttest{
#' dbs <- list(pgdb = LipidMS::pgdb, lysopgdb = LipidMS::lysopgdb,
#' fadb = LipidMS::fadb, adductsTable = LipidMS::adductsTable)
#'
#' candidates <- findCandidates(MS1 = LipidMS::mix_neg_fullMS, dbs[["pgdb"]],
#' ppm = 10, rt = c(min(MS1$RT), max(MS1$RT)), adducts = c("M-H"),
#' rttol = 3, dbs)
#'
#' MSMS <- rbind(LipidMS::mix_neg_Ce20, LipidMS::mix_neg_Ce40)
#' coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol)
#'
#' classConf <- checkClass(candidates, clfrags = c(227.0326, 209.022, 74.0359),
#' clrequired = c(F, F, F, F), ftype = c("F", "F", "NL"), ppm_products = 10,
#' dbs)
#'
#' sn1 <- chainFrags(coelfrags, chainfrags = c("lysopg_M-H"), ppm = 10,
#' dbs = dbs)
#' sn2 <- chainFrags(coelfrags, chainfrags = c("fa_M-H"), ppm = 10, dbs)
#'
#' chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
#'
#' intConf <- checkIntensityRules(intrules = c("lysopg_sn1/lysopg_sn1"),
#' rates = c("2/1"), intrequired = c(T), nchains=2, chainsComb, sn1, sn2)
#'
#' res <- organizeResults(candidates, clfrags = c(227.0326, 209.022, 74.0359),
#' classConf, chainsComb, intrules = c("lysopg_sn1/lysopg_sn1"), intConf,
#' nchains = 2, class="PG")
#' }
#'
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
organizeResults <- function(candidates, clfrags, classConf, chainComb, intrules,
                            intConf, nchains, class){
  results <- vector()
  if (nrow(candidates) > 0){
    if (nchains == 0){
      if (length(clfrags) > 0){
        idlevel <- "Subclass"
        confirmed <- classConf[[2]]
      } else {
        idlevel <- "MS-only"
        confirmed <- rep(TRUE, nrow(candidates))
      }
      results <- data.frame(ID = paste(class, "(", candidates$cb, ")", sep = ""),
                            Class = class,
                            CDB = candidates$cb, FAcomp = candidates$cb,
                            m.z = candidates$m.z, RT = candidates$RT,
                            I = candidates$int, Adducts = candidates$adducts,
                            ppm = candidates$ppms, confidenceLevel = idlevel,
                            stringsAsFactors = F)[confirmed,]
    } else if (nchains == 1){
      for (c in 1:nrow(candidates)){
        if (classConf[[2]][c] == T || length(clfrags) == 0){
          Class <- class
          CDB <- candidates$cb[c]
          m.z <- candidates$m.z[c]
          RT <- candidates$RT[c]
          I <- candidates$int[c]
          Adducts <- candidates$adducts[c]
          ppm <- candidates$ppm[c]
          ID <- paste(class, "(", CDB, ")", sep = "")
          FAcomp <- CDB
          if (chainComb[[c]] == T){
            confidenceLevel <- "FA"
            results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z, RT,
                                                 I, Adducts, ppm, confidenceLevel, stringsAsFactors = F))
          } else if (length(clfrags) > 0){
            confidenceLevel <- "Subclass"
            results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z, RT,
                                                 I, Adducts, ppm, confidenceLevel, stringsAsFactors = F))
          }
        }
      }
    } else if (nchains == 2){
      for (c in 1:nrow(candidates)){
        if (classConf[[2]][c] == T | length(clfrags) == 0){
          Class <- class
          CDB <- candidates$cb[c]
          m.z <- candidates$m.z[c]
          RT <- candidates$RT[c]
          I <- candidates$int[c]
          Adducts <- candidates$adducts[c]
          ppm <- candidates$ppm[c]
          if (nrow(chainComb[[c]]) > 0){
            for (comb in 1:nrow(chainComb[[c]])){
              if (intConf[[c]][comb] == T | length(intrules) == 0){
                ID <- paste(class, "(", paste(chainComb[[c]][comb,], collapse="/"),
                            ")", sep="")
                confidenceLevel <- "FA position"
              } else if (intConf[[c]][comb] == F | "Unknown" %in% intrules){
                ID <- paste(class, "(", paste(chainComb[[c]][comb,], collapse="_"),
                            ")", sep="")
                confidenceLevel <- "FA"
              }
              FAcomp <- paste(chainComb[[c]][comb,], collapse = " ")
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z,
                                                   RT, I, Adducts, ppm, confidenceLevel, stringsAsFactors = F))
            }
          } else {
            if (length(clfrags) > 0){
              ID <- paste(class, "(", CDB, ")", sep="")
              confidenceLevel <- "Subclass"
              FAcomp = ""
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z,
                                                   RT, I, Adducts, ppm, confidenceLevel, stringsAsFactors = F))
            }
          }
        }
      }
      if (!is.vector(results)){
        keep <- rep(NA, nrow(results))
        for (r in 1:nrow(results)){
          if (is.na(keep[r])){
            dup <- which(results$m.z == results$m.z[r] & results$RT ==
                           results$RT[r])
            dup <- dup[dup >= r]
            if (length(dup) > 1){
              rep <- duplicated(lapply(sapply(results$FAcomp[dup],
                                              strsplit, " "), sort))
              rep <- dup[which(rep)]
              all <- c(r, rep)
              tokeep <- all[which.max(LipidMS::confLevels[
                results$confidenceLevel[all], "order"])]
              if (results$confidenceLevel[tokeep] == "FA position"){
                tokeep <- all[which(results$confidenceLevel[all] ==
                                      results$confidenceLevel[tokeep])]
              }
              keep[tokeep] <- TRUE
              keep[setdiff(all, tokeep)] <- FALSE
            } else {
              keep[r] <- TRUE
            }
          }
        }
        results <- results[keep,]
        results <- unique(results)
      }
    } else if (nchains %in% c(3,4)){
      for (c in 1:nrow(candidates)){
        if (classConf[[2]][c] == T | length(clfrags) == 0){
          Class <- class
          CDB <- candidates$cb[c]
          m.z <- candidates$m.z[c]
          RT <- candidates$RT[c]
          I <- candidates$int[c]
          Adducts <- candidates$adducts[c]
          ppm <- candidates$ppm[c]
          if (nrow(chainComb[[c]]) > 0){
            for (comb in 1:nrow(chainComb[[c]])){
              if (intConf[[c]][comb] == T | length(intrules) == 0){
                ID <- paste(class, "(", paste(chainComb[[c]][comb,],
                                              collapse="/"), ")", sep="")
                confidenceLevel <- "FA position"
              } else if (intConf[[c]][comb] == F | intrules == "Unknown"){
                ID <- paste(class, "(", paste(chainComb[[c]][comb,],
                                              collapse="_"), ")", sep="")
                confidenceLevel <- "FA"
              }
              FAcomp <- paste(chainComb[[c]][comb,], collapse = " ")
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z, RT,
                                                   I, Adducts, ppm, confidenceLevel, stringsAsFactors = F))
            }
          } else {
            if (length(clfrags) > 0){
              ID <- paste(class, "(", CDB, ")", sep="")
              confidenceLevel <- "Subclass"
              FAcomp = ""
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z,
                                                   RT, I, Adducts, ppm, confidenceLevel, stringsAsFactors = F))
            }
          }
        }
      }
      if (!is.vector(results)){
        keep <- rep(NA, nrow(results))
        for (r in 1:nrow(results)){
          if (is.na(keep[r])){
            dup <- which(results$m.z == results$m.z[r] & results$RT ==
                           results$RT[r])
            dup <- dup[dup >= r]
            if (length(dup) > 1){
              rep <- duplicated(lapply(sapply(results$FAcomp[dup], strsplit,
                                              " "), sort))
              rep <- dup[which(rep)]
              all <- c(r, rep)
              tokeep <- all[which.max(LipidMS::confLevels[
                results$confidenceLevel[all], "order"])]
              if (results$confidenceLevel[tokeep] == "FA position"){
                tokeep <- all[which(results$confidenceLevel[all] ==
                                      results$confidenceLevel[tokeep])]
              }
              keep[tokeep] <- TRUE
              keep[setdiff(all, tokeep)] <- FALSE
            } else {
              keep[r] <- TRUE
            }
          }
        }
        results <- results[keep,]
        results <- unique(results)
      }
    }
  }
  return(results)
}
