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
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
combineChains <- function(candidates, nchains, sn1, sn2, sn3, sn4){
  if (nchains == 1){
    selected <- list()
    for (i in 1:nrow(candidates)){
      selected[[i]] <- candidates$cb[[i]] %in% sn1[[i]]$cb
    }
  }
  if (nchains == 2){
    chains1 <- list()
    chains2 <- list()
    for (c in 1:nrow(candidates)){
      if (nrow(sn1[[c]]) > 0){
        chains1[[c]] <- sn1[[c]]$cb
      } else {
        chains1[[c]] <- vector()
      }
      if (nrow(sn2[[c]]) > 0){
        chains2[[c]] <- sn2[[c]]$cb
      } else {
        chains2[[c]] <- vector()
      }
    }
    selected <- list()
    for (c in 1:nrow(candidates)){
      if(length(chains1[[c]]) > 0 & length(chains2[[c]]) > 0){
        x <- expand.grid(chains1[[c]], chains2[[c]], stringsAsFactors = F)
        comb <- list()
        for (i in 1:nrow(x)) {
          comb[[i]] <- paste(x[i, 1], x[i, 2], sep=" ")
        }
        sel <- unlist(lapply(comb, select, candidates$cb[c], n=2))
        selected[[c]] <- data.frame(sn1 = x[sel,1], sn2 = x[sel,2],
                                    stringsAsFactors = F)
      } else {
        selected[[c]] <- data.frame()
      }
    }
  }
  if (nchains == 3){
    chains1 <- list()
    chains2 <- list()
    chains3 <- list()
    for (c in 1:nrow(candidates)){
      if (nrow(sn1[[c]]) > 0){
        chains1[[c]] <- sn1[[c]]$cb
      } else {
        chains1[[c]] <- vector()
      }
      if (nrow(sn2[[c]]) > 0){
        chains2[[c]] <- sn2[[c]]$cb
      } else {
        chains2[[c]] <- vector()
      }
      if (nrow(sn3[[c]]) > 0){
        chains3[[c]] <- sn3[[c]]$cb
      } else {
        chains3[[c]] <- vector()
      }
    }
    selected <- list()
    for (c in 1:nrow(candidates)){
      if(length(chains1[[c]]) > 0 & length(chains2[[c]]) > 0 &
         length(chains3[[c]]) > 0){
        x <- expand.grid(chains1[[c]], chains2[[c]], chains3[[c]],
                         stringsAsFactors = F)
        comb <- list()
        for (i in 1:nrow(x)) {
          comb[[i]] <- paste(x[i, 1], x[i, 2], x[i, 3], sep=" ")
        }
        sel <- unlist(lapply(comb, select, candidates$cb[c], n=3))
        selected[[c]] <- data.frame(sn1 = x[sel,1], sn2 = x[sel,2],
                                    sn3 = x[sel,3], stringsAsFactors = F)
      } else {
        selected[[c]] <- data.frame()
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
        chains1[[c]] <- sn1[[c]]$cb
      } else {
        chains1[[c]] <- vector()
      }
      if (nrow(sn2[[c]]) > 0){
        chains2[[c]] <- sn2[[c]]$cb
      } else {
        chains2[[c]] <- vector()
      }
      if (nrow(sn3[[c]]) > 0){
        chains3[[c]] <- sn3[[c]]$cb
      } else {
        chains3[[c]] <- vector()
      }
      if (nrow(sn4[[c]]) > 0){
        chains4[[c]] <- sn4[[c]]$cb
      } else {
        chains4[[c]] <- vector()
      }
    }
    selected <- list()
    for (c in 1:nrow(candidates)){
      if(length(chains1[[c]]) > 0 & length(chains2[[c]]) > 0 &
         length(chains3[[c]]) > 0 & length(chains4[[c]]) > 0){
        x <- expand.grid(chains1[[c]], chains2[[c]], chains3[[c]], chains4[[c]],
                         stringsAsFactors = F)
        comb <- list()
        for (i in 1:nrow(x)) {
          comb[[i]] <- paste(x[i, 1], x[i, 2], x[i, 3], x[i, 4], sep=" ")
        }
        sel <- unlist(lapply(comb, select, candidates$cb[c], n=4))
        selected[[c]] <- data.frame(sn1 = x[sel,1], sn2 = x[sel,2],
                                    sn3 = x[sel,3], sn4 = x[sel,4],
                                    stringsAsFactors = F)
      } else {
        selected[[c]] <- data.frame()
      }
    }
  }
  return(selected)
}
