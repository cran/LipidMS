# combineChains
#' Combine chain fragments that could belong to the same precursor.
#'
#' It calculates combinations of chain fragments that sum up the same number of
#' carbons and double bounds as the precursor.
#'
#' @param candidates candidates data frame. Output of \link{findCandidates}.
#' @param nchains number of chains of the targeted lipid class.
#' @param sn1 list of chain fragments identified for sn1 position. Output of
#' \link{chainFrags}.
#' @param sn2 list of chain fragments identified for sn2 position. Output of
#' \link{chainFrags}. If required.
#' @param sn3 list of chain fragments identified for sn3 position. Output of
#' \link{chainFrags}. If required.
#' @param sn4 list of chain fragments identified for sn4 position. Output of
#' \link{chainFrags}. If required.
#'
#' @return List of data frames with candidate chains structures.
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
#' sn1 <- chainFrags(coelfrags, chainfrags = c("lysopg_M-H"), ppm = 10,
#' dbs = dbs)
#' sn2 <- chainFrags(coelfrags, chainfrags = c("fa_M-H"), ppm = 10, dbs = dbs)
#'
#' chainsComb <- combineChains(candidates, nchains=2, sn1, sn2)
#'
#' intConf <- checkIntensityRules(intrules = c("lysopg_sn1/lysopg_sn1"),
#' rates = c("2/1"), intrequired = c(T), nchains=2, chainsComb, sn1, sn2)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
combineChains <- function(candidates, nchains, sn1, sn2, sn3, sn4){
  if (nchains == 1){
    selected <- list()
    fragments <- rep(list(list()), nrow(candidates))
    for (c in 1:nrow(candidates)){
      fragments[[c]][[1]] <- sn1[[c]][which(sn1[[c]]$cb == candidates$cb[c]),]
      if (nrow(fragments[[c]][[1]]) > 0){
        selected[[c]] <- data.frame(sn1 = candidates$cb[[c]], stringsAsFactors = F)
      } else {
        selected[[c]] <- data.frame()
      }
    }
  }
  if (nchains == 2){
    chains1 <- list()
    chains2 <- list()
    for (c in 1:nrow(candidates)){
      if (nrow(sn1[[c]]) > 0){
        chains1[[c]] <- sn1[[c]]
      } else {
        chains1[[c]] <- data.frame()
      }
      if (nrow(sn2[[c]]) > 0){
        chains2[[c]] <- sn2[[c]]
      } else {
        chains2[[c]] <- data.frame()
      }
    }
    selected <- list()
    fragments <- rep(list(list()), nrow(candidates))
    for (c in 1:nrow(candidates)){
      if(length(chains1[[c]]) > 0 & length(chains2[[c]]) > 0){
        x <- expand.grid(unique(chains1[[c]]$cb), unique(chains2[[c]]$cb),
                         stringsAsFactors = F)
        comb <- list()
        for (i in 1:nrow(x)) {
          comb[[i]] <- paste(x[i, 1], x[i, 2], sep=" ")
        }
        sel <- unlist(lapply(comb, select, candidates$cb[c], n=2))
        selected[[c]] <- data.frame(sn1 = x[sel,1], sn2 = x[sel,2],
                                    stringsAsFactors = F)
        if (nrow(selected[[c]]) > 0){
          for (u in 1:nrow(selected[[c]])){
            fragments[[c]][[u]] <- unique(rbind(chains1[[c]][chains1[[c]]$cb %in%
                                                   selected[[c]]$sn1[u],],
                                    chains2[[c]][chains2[[c]]$cb %in%
                                                   selected[[c]]$sn2[u],]))
          }
        }
      } else {
        selected[[c]] <- data.frame()
        fragments[[c]] <- list()
      }
    }
  }
  if (nchains == 3){
    chains1 <- list()
    chains2 <- list()
    chains3 <- list()
    for (c in 1:nrow(candidates)){
      if (nrow(sn1[[c]]) > 0){
        chains1[[c]] <- sn1[[c]]
      } else {
        chains1[[c]] <- data.frame()
      }
      if (nrow(sn2[[c]]) > 0){
        chains2[[c]] <- sn2[[c]]
      } else {
        chains2[[c]] <- data.frame()
      }
      if (nrow(sn3[[c]]) > 0){
        chains3[[c]] <- sn3[[c]]
      } else {
        chains3[[c]] <- data.frame()
      }
    }
    selected <- list()
    fragments <- rep(list(list()), nrow(candidates))
    for (c in 1:nrow(candidates)){
      if(length(chains1[[c]]) > 0 & length(chains2[[c]]) > 0 &
         length(chains3[[c]]) > 0){
        x <- expand.grid(unique(chains1[[c]]$cb), unique(chains2[[c]]$cb),
                         unique(chains3[[c]]$cb), stringsAsFactors = F)
        comb <- list()
        for (i in 1:nrow(x)) {
          comb[[i]] <- paste(x[i, 1], x[i, 2], x[i, 3], sep=" ")
        }
        sel <- unlist(lapply(comb, select, candidates$cb[c], n=3))
        selected[[c]] <- data.frame(sn1 = x[sel,1], sn2 = x[sel,2],
                                    sn3 = x[sel,3], stringsAsFactors = F)
        if (nrow(selected[[c]]) > 0){
          for (u in 1:nrow(selected[[c]])){
            fragments[[c]][[u]] <- unique(rbind(chains1[[c]][chains1[[c]]$cb %in%
                                                        selected[[c]]$sn1[u],],
                                         chains2[[c]][chains2[[c]]$cb %in%
                                                        selected[[c]]$sn2[u],],
                                         chains3[[c]][chains3[[c]]$cb %in%
                                                        selected[[c]]$sn3[u],]))
          }
        }
      } else {
        selected[[c]] <- data.frame()
        fragments[[c]] <- list()
      }
    }
  }
  if (nchains == 4){
    chains1 <- list()
    chains2 <- list()
    chains3 <- list()
    chains4 <- list()
    for (c in 1:nrow(candidates)){
      if (nrow(sn1[[c]]) > 0){
        chains1[[c]] <- sn1[[c]]
      } else {
        chains1[[c]] <- vector()
      }
      if (nrow(sn2[[c]]) > 0){
        chains2[[c]] <- sn2[[c]]
      } else {
        chains2[[c]] <- vector()
      }
      if (nrow(sn3[[c]]) > 0){
        chains3[[c]] <- sn3[[c]]
      } else {
        chains3[[c]] <- vector()
      }
      if (nrow(sn4[[c]]) > 0){
        chains4[[c]] <- sn4[[c]]
      } else {
        chains4[[c]] <- vector()
      }
    }
    selected <- list()
    fragments <- rep(list(list()), nrow(candidates))
    for (c in 1:nrow(candidates)){
      if(length(chains1[[c]]) > 0 & length(chains2[[c]]) > 0 &
         length(chains3[[c]]) > 0 & length(chains4[[c]]) > 0){
        x <- expand.grid(unique(chains1[[c]]$cb), unique(chains2[[c]]$cb),
                         unique(chains3[[c]]$cb), unique(chains4[[c]]$cb),
                         stringsAsFactors = F)
        comb <- list()
        for (i in 1:nrow(x)) {
          comb[[i]] <- paste(x[i, 1], x[i, 2], x[i, 3], x[i, 4], sep=" ")
        }
        sel <- unlist(lapply(comb, select, candidates$cb[c], n=4))
        selected[[c]] <- data.frame(sn1 = x[sel,1], sn2 = x[sel,2],
                                    sn3 = x[sel,3], sn4 = x[sel,4],
                                    stringsAsFactors = F)
        if (nrow(selected[[c]]) > 0){
          for (u in 1:nrow(selected[[c]])){
            fragments[[c]][[u]] <- unique(rbind(chains1[[c]][chains1[[c]]$cb %in%
                                                        selected[[c]]$sn1[u],],
                                         chains2[[c]][chains2[[c]]$cb %in%
                                                        selected[[c]]$sn2[u],],
                                         chains3[[c]][chains3[[c]]$cb %in%
                                                        selected[[c]]$sn3[u],],
                                         chains4[[c]][chains4[[c]]$cb %in%
                                                        selected[[c]]$sn4[u],]))
          }
        }
      } else {
        selected[[c]] <- data.frame()
        fragments[[c]] <- list()
      }
    }
  }
  return(list(selected = selected, fragments = fragments))
}

