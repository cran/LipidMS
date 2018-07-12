# makexcmsSet
#' Convert peaklist from envipick to a xcmxSet to remove isotopes
#'
#' Convert peaklist from envipick to a xcmxSet to remove isotopes.
#'
#' @param envipick_object output of mzpick function from enviPick
#' @param files .mzXML file path
#' @param pD phenoData
#' @param polarity positive or negative
#'
#' @return xcmsSet
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
makexcmsSet <- function(envipick_object, files, pD, polarity){
  df <- envipick_object$Peaklist
  set1 <- as.matrix(cbind(mz = df[, "m/z"],
    mzmin = data.frame(df[, "m/z"]-df[, "var_m/z"]),
    mzmax = df[, "m/z"]+df[, "var_m/z"], rt = df[, "RT"],
    rtmin = df[, "minRT"], rtmax = df[, "maxRT"], into = df[, "sum_int"],
    intb = df[, "sum_int"]-df[, "sum_int"]/1000, maxo = df[, "max_int"],
    sn = rep(100, nrow(df)), sample = rep(1, nrow(df))))
  colnames(set1) <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into",
    "intb", "maxo", "sn", "sample")
  object <- methods::new("xcmsSet")
  object@peaks <- set1
  object@phenoData <- data.frame(pD, row.names = files)
  colnames(object@phenoData) <- "class"
  object@rt <- list(raw = envipick_object$Scans[1],
    corrected = envipick_object$Scans[1])
  object@filepaths <- paste(c(getwd(), files), collapse = "/")
  object@profinfo <- list(method = "bin", step = 0.1)
  object@polarity <- polarity

  return(object)
}
