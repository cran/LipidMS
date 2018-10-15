# organizeResults
#' Prepare output for LipidMS annotation functions
#'
#' Prepare a readable output for LipidMS identification functions.
#'
#' @param candidates candidates data frame. Output of \link{findCandidates}.
#' @param clfrags vector containing the expected fragments for a given lipid
#' class.
#' @param classConf output of \link{checkClass}
#' @param chainsComb output of \link{combineChains}
#' @param intrules character vector specifying the fragments to compare. See
#' \link{checkIntensityRules}.
#' @param intConf output of \link{checkIntensityRules}
#' @param nchains number of chains of the targeted lipid class.
#' @param class character value. Lipid class (i.e. PC, PE, DG, TG, etc.).
#'
#' @examples
#' \donttest{
#' dbs <- assignDB()
#'
#' candidates <- findCandidates(MS1 = LipidMS::MS1_neg$peaklist,
#' db = dbs$pgdb, ppm = 10, rt = c(0, 2000), adducts = c("M-H"),
#' rttol = 10, rawData = MS1_neg$rawScans, coelCutoff = 0.8)
#'
#' MSMS <- rbind(LipidMS::MSMS1_neg$peaklist, LipidMS::MSMS2_neg$peaklist)
#' rawData <- rbind(LipidMS::MS1_neg$rawScans, LipidMS::MSMS1_neg$rawScans,
#' LipidMS::MSMS2_neg$rawScans)
#' coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol = 10, rawData = rawData,
#' coelCutoff = 0.8)
#'
#' classConf <- checkClass(candidates, coelfrags,
#' clfrags = c(227.0326, 209.022, 74.0359), clrequisites = c(F, F, F, F),
#' ftype = c("F", "F", "NL"), ppm = 10, dbs = dbs)
#'
#' sn1 <- chainFrags(coelfrags, chainfrags = c("lysopg_M-H"), ppm = 10,
#' dbs = dbs)
#' sn2 <- chainFrags(coelfrags, chainfrags = c("fa_M-H"), ppm = 10, dbs = dbs)
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
organizeResults <- function(candidates, clfrags, classConf, chainsComb, intrules,
                            intConf, nchains, class){
  results <- vector()
  if (nrow(candidates) > 0){
    if (nchains == 0){
      for (c in 1:nrow(candidates)){
        if (length(clfrags) > 0){
          idlevel <- "Subclass"
          if (classConf$passed[c] == T){
            score <- mean(classConf$fragments[[c]]$coelScore, na.rm = T)
          }
        } else {
          idlevel <- "MS-only"
          score <- 0
        }
        if (classConf$passed[c] == T ||
            classConf$passed[c] == "No class fragments"){
          results <- rbind(results,
                           data.frame(ID = paste(class, "(", candidates$cb[c], ")",
                                                 sep = ""),Class = class,
                                      CDB = candidates$cb[c],
                                      FAcomp = candidates$cb[c],
                                      m.z = candidates$m.z[c],
                                      RT = candidates$RT[c],
                                      I = candidates$int[c],
                                      Adducts = candidates$adducts[c],
                                      ppm = round(candidates$ppms[c], 3),
                                      confidenceLevel = idlevel,
                                      peakID = candidates$peakID[c],
                                      PFCS = round(score, 3),
                                      stringsAsFactors = F))
        }
      }
    } else if (nchains == 1){
      for (c in 1:nrow(candidates)){
        if (classConf$passed[c] == T || length(clfrags) == 0){
          Class <- class
          CDB <- candidates$cb[c]
          m.z <- candidates$m.z[c]
          RT <- candidates$RT[c]
          I <- candidates$int[c]
          Adducts <- candidates$adducts[c]
          ppm <- round(candidates$ppm[c], 3)
          ID <- paste(class, "(", CDB, ")", sep = "")
          FAcomp <- CDB
          if (nrow(chainsComb$selected[[c]]) > 0){
            confidenceLevel <- "FA"
            if (length(clfrags) > 0){
              score <-  mean(c(classConf$fragments[[c]]$coelScore,
                            chainsComb$fragments[[c]]$coelScore, na.rm = T))
            } else {
              score <-  mean(chainsComb$fragments[[c]]$coelScore, na.rm = T)
            }
            results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z, RT,
                                                 I, Adducts, ppm, confidenceLevel,
                                                 peakID = candidates$peakID[c],
                                                 PFCS = round(score, 3),
                                                 stringsAsFactors = F))
          } else if (length(clfrags) > 0){
            confidenceLevel <- "Subclass"
            score <- mean(classConf$fragments[[c]]$coelScore, na.rm = T)
            results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z, RT,
                                                 I, Adducts, ppm, confidenceLevel,
                                                 peakID = candidates$peakID[c],
                                                 PFCS = round(score, 3),
                                                 stringsAsFactors = F))
          }
        }
      }
    } else if (nchains == 2){
      for (c in 1:nrow(candidates)){
        if (classConf$passed[c] == T || length(clfrags) == 0){
          Class <- class
          CDB <- candidates$cb[c]
          m.z <- candidates$m.z[c]
          RT <- candidates$RT[c]
          I <- candidates$int[c]
          Adducts <- candidates$adducts[c]
          ppm <- round(candidates$ppm[c], 3)
          if (nrow(chainsComb$selected[[c]]) > 0){
            for (comb in 1:nrow(chainsComb$selected[[c]])){
              if (intConf[[c]][comb] == T || length(intrules) == 0){
                ID <- paste(class, "(", paste(chainsComb$selected[[c]][comb,],
                                              collapse="/"),
                            ")", sep="")
                confidenceLevel <- "FA position"
              } else if (intConf[[c]][comb] == F || "Unknown" %in% intrules){
                ID <- paste(class, "(", paste(chainsComb$selected[[c]][comb,],
                                              collapse="_"),
                            ")", sep="")
                confidenceLevel <- "FA"
              }
              if (length(clfrags) > 0){
                score <-  mean(c(classConf$fragments[[c]]$coelScore,
                                 chainsComb$fragments[[c]][[comb]]$coelScore),
                               na.rm = T)
              } else {
                score <-  mean(chainsComb$fragments[[c]][[comb]]$coelScore, na.rm = T)
              }
              FAcomp <- paste(chainsComb$selected[[c]][comb,], collapse = " ")
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z,
                                                   RT, I, Adducts, ppm,
                                                   confidenceLevel,
                                                   peakID = candidates$peakID[c],
                                                   PFCS = round(score, 3),
                                                   stringsAsFactors = F))
            }
          } else {
            if (length(clfrags) > 0){
              ID <- paste(class, "(", CDB, ")", sep="")
              confidenceLevel <- "Subclass"
              FAcomp = ""
              score <- mean(classConf$fragments[[c]]$coelScore, na.rm = T)
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z,
                                                   RT, I, Adducts, ppm,
                                                   confidenceLevel,
                                                   peakID = candidates$peakID[c],
                                                   PFCS = round(score, 3),
                                                   stringsAsFactors = F))
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
        if (classConf$passed[c] == T || length(clfrags) == 0){
          Class <- class
          CDB <- candidates$cb[c]
          m.z <- candidates$m.z[c]
          RT <- candidates$RT[c]
          I <- candidates$int[c]
          Adducts <- candidates$adducts[c]
          ppm <- round(candidates$ppm[c], 3)
          if (nrow(chainsComb$selected[[c]]) > 0){
            for (comb in 1:nrow(chainsComb$selected[[c]])){
              if (intConf[[c]][comb] == T || length(intrules) == 0){
                ID <- paste(class, "(", paste(chainsComb$selected[[c]][comb,],
                                              collapse="/"), ")", sep="")
                confidenceLevel <- "FA position"
              } else if (intConf[[c]][comb] == F || intrules == "Unknown"){
                ID <- paste(class, "(", paste(chainsComb$selected[[c]][comb,],
                                              collapse="_"), ")", sep="")
                confidenceLevel <- "FA"
              }
              if (length(clfrags) > 0){
                score <-  mean(c(classConf$fragments[[c]]$coelScore,
                                 chainsComb$fragments[[c]][[comb]]$coelScore),
                               na.rm = T)
              } else {
                score <-  mean(chainsComb$fragments[[c]][[comb]]$coelScore,
                               na.rm = T)
              }
              FAcomp <- paste(chainsComb$selected[[c]][comb,], collapse = " ")
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z, RT,
                                                   I, Adducts, ppm,
                                                   confidenceLevel,
                                                   peakID = candidates$peakID[c],
                                                   PFCS = round(score, 3),
                                                   stringsAsFactors = F))
            }
          } else {
            if (length(clfrags) > 0){
              ID <- paste(class, "(", CDB, ")", sep="")
              confidenceLevel <- "Subclass"
              FAcomp = ""
              score <- mean(classConf$fragments[[c]]$coelScore, na.rm = T)
              results <- rbind(results, data.frame(ID, Class, CDB, FAcomp, m.z,
                                                   RT, I, Adducts, ppm,
                                                   confidenceLevel,
                                                   peakID = candidates$peakID[c],
                                                   PFCS = round(score, 3),
                                                   stringsAsFactors = F))
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
