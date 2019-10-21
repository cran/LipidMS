# checkClass
#' Search of class fragments to confirm the lipid class.
#'
#' Search of characteristic fragments that confirm a given lipid class.
#'
#' @param candidates output of \link{findCandidates} function.
#' @param coelfrags list of peaks coeluting with each candidate. Output of
#' \link{coelutingFrags}.
#' @param clfrags vector containing the expected fragments for a given lipid
#' class. See details.
#' @param ftype character vector indicating the type of fragments in clfrags.
#' It can be: "F" (fragment), "NL" (neutral loss) or "BB" (building block).
#' See details.
#' @param clrequisites logical vector indicating if each class fragment is
#' required or not. If none of the fragment is required, at least one of them
#' must be present within the coeluting fragments. If the presence of any
#' fragment excludes the class, it can be specified by using "excluding".
#' @param ppm m/z tolerance in ppm.
#' @param dbs list of data bases required for the annotation. By default, dbs
#' contains the required data frames based on the default fragmentation rules.
#' If these rules are modified, dbs may need to be changed. If data bases have
#' been customized using \link{createLipidDB}, they also have to be modified
#' here. It is employed when some fragment belongs to "BB" \code{ftype}.
#'
#'
#' @details \code{clfrags}, \code{ftype} and \code{clrequisites} will indicate
#' the rules to confirm a lipid class. All three arguments must have the same
#' length.
#'
#' This function allows three different types of fragments: fragments with a
#' specific m/z as for example 227.0326 for PG in negative mode, which needs to
#' be defined as clfrags = c(227.0326) and ftype = c("F"); neutral losses such
#' as the head group of some PL (i.e. NL of 74.0359 in PG in negative mode),
#' which will be defined as clfrags = c(74.0359) and ftype = c("NL"); or
#' building blocks resulting from the loss of some groups, as for example, PA as
#' M-H resulting from the loss of the head group (glycerol) in PG in ESI-, which
#' will be defined as clfrags = c("pa_M-H") and ftype = c("BB"). The last
#' two options could define the same fragments. In this case just one of them
#' would be necessary.
#'
#' When using the third type of fragment ("BB"), the building block will be
#' specified in lower case (i.e. pa, dg, lysopa, mg, etc.) and the adduct will be
#' given as it appears in the adductsTable, both separated by "_". Names for the
#' building blocks are the ones used for the LipidMS databases without the "db"
#' at the end.
#'
#' In case the presence of a fragment indicates that the candidate does not
#' belong to the lipid class (i.e. loss of CH3 in PE, which corresponds to a PC
#' actually), this will be specified by using \code{clrequisites} = c("excluding").
#'
#' @return List with 2 elements: a matrix with logical values (presence/absense)
#' of each expected fragment (columns) for each candidate (rows), and a logical
#' vector with the confirmation of the lipid class for each candidate.
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
#' rawData <- rbind(MS1_neg$rawScans, MSMS1_neg$rawScans,
#' MSMS2_neg$rawScans)
#' coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol = 10, rawData = rawData,
#' coelCutoff = 0.8)
#'
#' classConf <- checkClass(candidates, coelfrags,
#' clfrags = c(227.0326, 209.022, 74.0359), clrequisites = c(F, F, F, F),
#' ftype = c("F", "F", "NL"), ppm = 10, dbs = dbs)
#' }
#' \donttest{
#' library(LipidMSdata)
#' dbs <- assignDB()
#'
#' candidates <- findCandidates(MS1 = MS1_neg$peaktable,
#' db = dbs$pgdb,
#' ppm = 10, rt = c(0, 2000), adducts = c("M-H"),
#' rttol = 10, rawData = MS1$rawData, coelCutoff = 0.8)
#'
#' MSMS <- rbind(MSMS1_neg$peaktable, MSMS2_neg$peaktable)
#' rawData <- rbind(MS1_neg$rawData, MSMS1_neg$rawData,
#' MSMS2_neg$rawData)
#' coelfrags <- coelutingFrags(candidates$RT, MSMS, rttol = 10, rawData = rawData,
#' coelCutoff = 0.8)
#'
#' classConf <- checkClass(candidates, coelfrags,
#' clfrags = c(227.0326, 209.022, 74.0359), clrequisites = c(F, F, F, F),
#' ftype = c("F", "F", "NL"), ppm = 10, dbs = dbs)
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
checkClass <- function(candidates, coelfrags, clfrags, ftype, clrequisites,
                       ppm = 10, dbs){
  adductsTable <- dbs$adductsTable
  allfragments <- rep(list(data.frame()), nrow(candidates))
  if (length(clfrags) == 0){
    return(list(presence = matrix(), passed = rep("No class fragments",
                                                  nrow(candidates))))
  } else {
    verified <- matrix(NA, ncol=length(clfrags), nrow = nrow(candidates))
    for (i in 1:length(clfrags)){
      if (ftype[i] == "F"){
        classf <- lapply(coelfrags, function(x){
          filtermsms(x, as.numeric(clfrags[i]), ppm)})
        verified[,i] <- unlist(lapply(classf, function(x) {if(is.data.frame(x)){nrow(x) > 0}else {FALSE}}))
        allfragments <- Map(rbind, classf, allfragments)
      } else if (ftype[i] == "NL") {
        classf <- list()
        for (c in 1:nrow(candidates)){
          f <- filtermsms(coelfrags[[c]], candidates$m.z[c]-as.numeric(clfrags[i]), ppm)
          classf[[c]] <- f
        }
        verified[,i] <- unlist(lapply(classf, function(x) {if(is.data.frame(x)){nrow(x) > 0}else {FALSE}}))
        allfragments <- Map(rbind, classf, allfragments)
      } else if (ftype[i] == "BB"){
        db <- dbs[[paste(unlist(strsplit(clfrags[i], "_"))[1], "db", sep="")]]
        mdiff <- adductsTable[adductsTable$adduct ==
                                unlist(strsplit(clfrags[i], "_"))[2], "mdiff"]
        classf <- list()
        for (c in 1:nrow(candidates)){
          f <- filtermsms(coelfrags[[c]],
                          db$Mass[which(db$total == candidates$cb[c])]+ mdiff,
                          ppm)
          classf[[c]] <- f
        }
        verified[,i] <- unlist(lapply(classf, function(x) {if(is.data.frame(x)){nrow(x) > 0}else {FALSE}}))
        allfragments <- Map(rbind, classf, allfragments)
      }
    }
    if (any(clrequisites == T | clrequisites == "TRUE")){ # when there are required fragments, we check them
      passed <- apply(verified, 1,
                      function(x) sum(x == T & (clrequisites == T |
                                                  clrequisites == "TRUE"))) ==
        sum(clrequisites == T | clrequisites == "TRUE")
    } else { # if there isnt any required fragment, any of them will be enough
      passed <- apply(verified, 1, sum) > 0
    }
    if (any(clrequisites == "excluding")){
      excl <- which(clrequisites == "excluding")
      for (e in 1:length(excl)){
        presence = verified[,excl[e]]
        passed[presence == T] = FALSE
      }
    }
    return(list(presence = verified, passed = passed, fragments = allfragments))
  }
}





