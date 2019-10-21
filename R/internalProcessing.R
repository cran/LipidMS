# annotateIsotopes
#' Annotate isotopes
#'
#' Annotate isotopes based on mass differences, retention time and peak
#' correlation if required.
#'
#' @param peaklist extracted peaks. Data.frame with 4 columns (m.z, RT, int
#' and peakID).
#' @param rawScans raw scan data. Data.frame with 5 columns (m.z, RT, int,
#' peakID and Scan).
#' @param ppm mass tolerance in ppm.
#' @param rttol rt windows with the same units used in peaklist.
#'
#' @return peaklist with 6 columns (m.z, RT, int, peakID, isotope and group).
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
annotateIsotopes <- function(peaklist, rawScans, ppm, rttol){
  peaklist <- peaklist[order(peaklist[,"m.z"], decreasing = F),]
  anno_peaklist <- peaklist[c(),]
  cluster <- 1 # contador
  while (nrow(peaklist) > 0){ # iniciamos un bucle para pasar los clusters agrupados
    # a una nueva peaklist, hasta que la anterior se vacie
    a <- peaklist[1,]
    colnames(a) <- c("m.z", "RT", "int", "peakID")
    ss <- peaklist[which(abs(peaklist[,"RT"] - a[,"RT"]) < rttol &
                           peaklist[,"m.z"] >= a[,"m.z"]),]
    ss <- ss[order(ss[,"m.z"]),]
    if (nrow(ss) > 1){
      dm <- (ss[,"m.z"]-a[,"m.z"])
      ss[,"iso"] <- round(dm,0)
      errors <- abs(dm-ss[,"iso"]*1.003355)*1e6/(a[,"m.z"]+round(dm, 0)*1.003355)
      group <- cbind(ss[errors < ppm,c("m.z", "RT", "int", "peakID", "iso"),])
      if (nrow(group) > 1){
        keep <- c(TRUE, rep(FALSE, nrow(group)-1))
        for (i in 2:nrow(group)){
          prev <- which(group[,"iso"] - group[i,"iso"] == -1 & keep[i-1] == TRUE)
          if(length(prev) > 0 & group[i,"iso"] == 1){
            intcheck <- any((group[i,"int"]/group[prev, "int"]) > 0.011 &
                              (group[i,"int"]/group[prev, "int"]) <
                              group[prev[1],"m.z"]*0.011/12)
            keep[i] <- intcheck
          } else if(length(prev) > 0 & group[i,"iso"] > 1){
            keep[i] <- any(group[,"iso"] - group[i,"iso"] == -1 & keep[i-1] == TRUE)
          } else {
            break
          }
        }
        group <- group[keep,]
        keep <- c(TRUE, rep(TRUE, nrow(group)-1))
        check_overlaps <- which(table(group[,"iso"]) > 1)
        if (length(check_overlaps) > 0){
          repeated <- as.numeric(names(check_overlaps))
          cors <- sapply(group[group[,"iso"] %in% repeated,"peakID"],
                         coelutionScore, group[1,"peakID"], rawScans)
          remove <- which(cors < 0.8)
          if (length(remove) > 0 & length(remove) < length(repeated)){
            keep[which(group[,"iso"] %in% repeated)[remove]] <- FALSE
          } else if (length(remove) > 0 & length(remove) == length(repeated)){
            remove[which.max(cors)] <- FALSE
            keep[which(group[,"iso"] %in% repeated)[remove]] <- FALSE
          }
        }
        group <- group[keep,]
      }
      if ( nrow(group) > 1){
        group[,"cluster"] <- cluster
        cluster <- cluster + 1
      } else {
        group[,"cluster"] <- 0
      }
    } else {
      ss[,"iso"] <- 0
      group <- ss
      group[,"cluster"] <- 0
    }
    anno_peaklist <- rbind(anno_peaklist, group)
    peaklist <- peaklist[which(!peaklist[,"peakID"] %in% group[,"peakID"]),]
  }
  anno_peaklist[,"iso"] <- paste("[M+", anno_peaklist[,"iso"], "]", sep="")
  anno_peaklist[,"iso"][anno_peaklist[,"cluster"] == 0] <- ""
  colnames(anno_peaklist) <- c("m.z", "RT", "int", "peakID", "isotope", "group")
  return(anno_peaklist)
}


