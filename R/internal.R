# sumChains
#' Calculate total number of carbons and double bounds of lipid chains
#'
#' Given the structure of a lipid specie, it sums up the chains.
#'
#' @param chains character value with the configuration of the chains separated
#' by a white space
#' @param n number of chains
#'
#' @return Character value indicating the total number of carbons and double
#' bounds
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
sumChains <- function(chains, n){
  cb <- unlist(strsplit(chains, "[: ]"))
  if (n == 2){
    sum_c <- as.numeric(cb[1])+as.numeric(cb[3])
    sum_b <- as.numeric(cb[2])+as.numeric(cb[4])
  } else if (n == 3){
    sum_c <- as.numeric(cb[1])+as.numeric(cb[3])+as.numeric(cb[5])
    sum_b <- as.numeric(cb[2])+as.numeric(cb[4])+as.numeric(cb[6])
  } else if (n == 4){
    sum_c <- as.numeric(cb[1])+as.numeric(cb[3])+as.numeric(cb[5])+as.numeric(cb[7])
    sum_b <- as.numeric(cb[2])+as.numeric(cb[4])+as.numeric(cb[6])+as.numeric(cb[8])
  }
  total <- paste(c(sum_c, sum_b), collapse=":")
  return(total)
}

# mzMatch
#' mz match withing a vector of mz values
#'
#' This function searches marches between a given mz and a vector of mz values
#' with certain mass  tolerance and returns the index of the matched values. It
#' is used by identification functions to find candidates of each class of lipid
#' based on full MS information.
#'
#' @param mz mz value to be matched
#' @param mzvector vector of mz values
#' @param ppm mass error tolerance
#'
#' @return Numeric vector indicating the index of matched mz values and ppms for
#' each one of those matches (match1, ppm1, match2, ppm2, etc.)
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
mzMatch <- function(mz, mzvector, ppm){
  matches_ppm <- vector()
  for (i in 1:length(mzvector)){
    ppm_observed <- abs(((mz-mzvector[i])/mz)*1000000)
    if (ppm_observed <= ppm){
      matches_ppm <- append(matches_ppm, c(i, ppm_observed))
    }
  }
  return(matches_ppm)
}

# cbs
#' Total number of carbons and double bounds
#'
#' This function matches mz values with neutral masses from a dataframe which
#' links masses and structures (carbons and double bounds) and extracts the
#' structural information. It is used by identification functions to look for
#' the structure of the previously chosen candidates.
#'
#' @param mz mz value to be matched
#' @param ppm mass error tolerance
#' @param db database
#' @param charge numeric value indicating the charge of the ion
#'
#' @return Character value or vector indicating structural information
#' (carbons:bounds)
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
cbs <- function(mz, ppm, db, charge=0){
  cb <- as.vector(db[abs(db["Mass"]+charge-mz) <
      ppm*(db["Mass"]+charge)/1000000, "total"])
  if (length(cb) > 1){
    cb <- cb[1]
  }
  if (length(cb) == 0){
    cb <- ""
  }
  return(cb)
}

# filtermsms
#' Presence or absence of an mz value withing a vector of mz values
#'
#' This function indicates the presence or absence of a fragment within a set
#' of mz values with certain tolerance. It is used by identification functions
#' to look for the generic fragments of each class of lipid.
#'
#' @param fragments vector of mz values
#' @param frag mz to be matched
#' @param ppm mass tolerance
#'
#' @return Logical value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
filtermsms <- function(fragments, frag, ppm){
  if (sum(abs((fragments - frag)/frag)*1000000 < ppm) > 0){
    res <- TRUE
  } else {
    res <- FALSE
  }
  return(res)
}

# frags
#' Search for fragments of interest withing a list of coeluting fragments
#'
#' Given a set of coeluting fragments, this function searches for matches within
#' a database. It is used by identification functions to extract fragments of
#' interest based on the fragmentation patterns of each class of lipid.
#'
#' @param df data frame containing coeluting fragments
#' @param ppm mass tolerance
#' @param db database (data frame with two columns) where to look into
#' @param charge mdiff
#'
#' @return Data frame containing matched ions information
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
frags <- function(df, ppm, db, charge){
  if (nrow(df) > 0){
    cb <- data.frame(0,0,0,0)
    colnames(cb) <- c("cb", "m.z", "RT", "int")
    found <- FALSE
    if (nrow(df) > 0){
      for (x in 1:nrow(df)){
        y <- abs((db[,"Mass"]+charge-df[x,"m.z"])*1000000/
                   (db[,"Mass"]+charge)) < ppm
        if (sum(y) > 0){
          cb <- rbind(cb, c(db[which(y == TRUE), "total"],
                            df[x,"m.z"], df[x,"RT"], df[x,"int"]))
          found <- TRUE
        }
      }
      if (found == FALSE){
        cb <- data.frame()
      } else {
        cb <- cb[2:nrow(cb),]
        if (sum(duplicated(cb$cb)) > 0){
          cb <- joinfrags(cb)
        }
      }
    }
  } else {
    cb <- data.frame()
  }
  return(cb)
}

# joinfrags
#' Join fragments information when several peaks of the same fragment are
#' coeluting with a unique candidate
#'
#' Function employed by \link{frags}.
#'
#' @param df data frame containing coeluting fragments
#'
#' @return Data frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
joinfrags <- function(df){
  new <- vector()
  for (f in unique(df$cb)){
    new <- rbind(new, data.frame(cb=f, m.z=mean(as.numeric(df$m.z[df$cb == f])),
                                 RT = mean(as.numeric(df$RT[df$cb == f])),
                                 int = sum(as.numeric(df$int[df$cb == f])),
                                 stringsAsFactors = F))
  }
  return(new)
}

# diffcb
#' Difference between two carbon:bounds structures
#'
#' This function calculates the number of carbon and double bounds that differ
#' between two structures.
#'
#' @param total character value indicating the precursor structure
#' @param frag character value indicating the fragment structure
#' @param db db of chains to be considered
#'
#' @return Character value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
diffcb <- function(total, frag, db){
  fas <- db$total
  total <- unlist(strsplit(total, ":"))
  fr <- unlist(sapply(frag, strsplit, ":"))
  x <- as.numeric(total) - as.numeric(fr)
  frags <- paste(x[seq(1, length(x), 2)], x[seq(2, length(x), 2)], sep=":")
  del <- setdiff(frags, fas)
  frags[frags %in% del] <- ""
  if (length(frags) > 0){
    return(as.vector(frags))
  } else {
    return("")
  }
}

# select
#' Check matches between chains composition and precursor structures
#'
#' This function checks if the sum up of the chains structure match the
#' precursor structure. It is used by \link{combineChains}.
#'
#' @param chains character value containing chains structure separated by a
#' white space.
#' @param parent precursor ion structure
#' @param n number of chains
#'
#' @return Logical value
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
select <- function (chains, parent, n){
  chains <- unlist(strsplit(chains, "[: ]"))
  if (n == 3) {
    sum_c <- as.numeric(chains[1]) + as.numeric(chains[3]) +
      as.numeric(chains[5])
    sum_e <- as.numeric(chains[2]) + as.numeric(chains[4]) +
      as.numeric(chains[6])
  } else if (n == 2) {
    sum_c <- as.numeric(chains[1]) + as.numeric(chains[3])
    sum_e <- as.numeric(chains[2]) + as.numeric(chains[4])
  } else if (n == 4){
    sum_c <- as.numeric(chains[1])+as.numeric(chains[3])+as.numeric(chains[5])+
      as.numeric(chains[7])
    sum_e <- as.numeric(chains[2])+as.numeric(chains[4])+as.numeric(chains[6])+
      as.numeric(chains[8])
  }
  equal <- paste(c(sum_c, sum_e), collapse = ":") == parent
  return(equal)
}

# findPrecursor
#' Find candidate precursor from fullMS function
#'
#' This function is employed by all identification function in this package to
#' find possible precursors in the fullMS function.
#'
#' @param MS1 data frame containing m/z, RT and intensity of ions from the first
#' function
#' @param db database to be searched for matches
#' @param massdif mass difference between neutral mass and the adduct expected
#' @param rt rt window
#' @param n numeric value indicating whether to look for monomers (1), dimers
#' (2), etc.
#' @param charge numeric value indicating the expected charge of the ions
#'
#'
#' @return Subset of the original data frame without adducts
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
findPrecursor <- function(MS1, db, ppm, massdif, rt, n=1, charge=1){
  precursors <- unlist(lapply((n*db$Mass+massdif)/abs(charge), mzMatch,
                              MS1$m.z, ppm))
  if (length(precursors) > 0){
    matches <- precursors[seq(1, length(precursors), 2)]
    ppms <- precursors[seq(2, length(precursors), 2)]
    prec <- MS1[matches,]
    if (class(prec) == "numeric"){
      prec <- as.data.frame(t(prec))
    }
    ppms <- ppms[prec$RT >= rt[1] & prec$RT <= rt[2]]
    prec <- prec[prec$RT >= rt[1] & prec$RT <= rt[2],]
    if (nrow(prec) > 0){
      if (n == 1){
        cb <- sapply(prec$m.z*abs(charge), cbs, ppm+2, db, massdif)
      } else if (n == 2){
        cb <- vector()
        for (i in 1:nrow(prec)){
          cb <- append(cb, cbs((prec$m.z[i]*abs(charge)-massdif)/2, ppm+2, db))
        }
      }
      # joining info
      allprec <- cbind(prec, ppms, as.data.frame(cb, stringsAsFactors = F))
      return(allprec)
    } else {
      return(data.frame())
    }
  } else {
    return(data.frame())
  }
}


# joinAdducts
#' Join information about precursor candidates found using different adducts
#'
#' This function joins tables of precursor candidates found using different
#' adducts.
#'
#' @param df1 data frame containing identification results using the main adduct
#' @param df data frame containing identification results using a secondary
#' adduct
#' @param rttol retention time tolerance in seconds
#' @param main character value indicating the preferred adduct
#' @param other character value indicating the ceecondary adduct
#'
#' @return Subset of the original data frame without adducts
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
joinAdducts <- function(df1, df2, rttol, main, other){
  if (nrow(df2) > 0 & nrow(df1) == 0){
    df1 <- df2
    df2 <- data.frame()
    main <- other
  }
  if (nrow(df1) > 0 & nrow(df2) > 0){
    mzmatchesdf1 <- vector()
    mzmatchesdf2 <- vector()
    for (i in 1:nrow(df2)){
      m <- which(abs(df1$m.z - df2$m.z[i]) < 0.001 & abs(df1$RT - df2$RT[i]) <
                   rttol/2)
      mzmatchesdf1 <- append(mzmatchesdf1, m)
      if (length(m) > 0){
        mzmatchesdf2 <- append(mzmatchesdf2, i)
      }
    }
    df3 <- df2[mzmatchesdf2,]
    if (length(mzmatchesdf1) > 0){
      df4 <- df2[-mzmatchesdf2,]
    } else {
      df4 <- df2
    }
    if (nrow(df4) > 0){
      cbmatchesdf1 <- vector()
      cbmatchesdf4 <- vector()
      for (i in 1:nrow(df4)){
        m <- which(df1$cb == df4$cb[i] & abs(df1$RT - df4$RT[i]) < rttol)
        cbmatchesdf1 <- append(cbmatchesdf1, m)
        if (length(m) > 0){
          cbmatchesdf4 <- append(cbmatchesdf4, i)
        }
      }
      if ("adducts" %in% colnames(df1)){
        adducts <- as.vector(df1$adducts)
        adducts[cbmatchesdf1] <- paste(adducts[cbmatchesdf1], other, sep = ";")
      } else {
        adducts <- rep(main, nrow(df1))
        adducts[cbmatchesdf1] <- paste(main, other, sep = ";")
      }
      df1$adducts <- adducts
      dif <- setdiff(c(1:nrow(df4)), cbmatchesdf4)
      if (length(dif) > 0){
        df1 <- rbind(df1, cbind(df4[dif,], adducts=rep(other, length(dif))))
      }
    }
    if (nrow(df3) > 0){
      cbmatchesdf1 <- vector()
      cbmatchesdf3 <- vector()
      for (i in 1:nrow(df3)){
        m <- which(df1$cb == df3$cb[i] & abs(df1$RT - df3$RT[i]) < rttol)
        cbmatchesdf1 <- append(cbmatchesdf1, m)
        if (length(m) > 0){
          cbmatchesdf3 <- append(cbmatchesdf3, i)
        }
      }
      if ("adducts" %in% colnames(df1)){
        adducts <- as.vector(df1$adducts)
        adducts[cbmatchesdf1] <- paste(adducts[cbmatchesdf1], other, sep = ";")
      } else {
        adducts <- rep(main, nrow(df1))
        adducts[cbmatchesdf1] <- paste(main, other, sep = ";")
      }
      df1$adducts <- adducts
    }
  } else {
    if (!"adducts" %in% colnames(df1)){
      df1$adducts <- rep(main, nrow(df1))
    }
  }
  return(df1)
}
# checkIntRules
#' Check intensity rules
#'
#' Check intensity rules to confirm chains structure.
#'
#' @param intrules character vector specifying the fragments to compare. See
#' details.
#' @param rates character vector with the expected rates given as a string
#' (i.e. "3/1"). See details.
#' @param intrequired logical vector indicating if any of the rules is required.
#' If not, at least one must be verified to confirm the structure.
#' @param nchains number of chains of the targeted lipid class.
#' @param combinations output of \link{combineChains}
#' @param sn1 list of chain fragments identified for sn1 position. Output of
#' \link{chainFrags}.
#' @param sn2 list of chain fragments identified for sn2 position. Output of
#' \link{chainFrags}. If required.
#' @param sn3 list of chain fragments identified for sn3 position. Output of
#' \link{chainFrags}. If required.
#' @param sn4 list of chain fragments identified for sn4 position. Output of
#' \link{chainFrags}. If required.
#'
#' @details This function will be employed when the targeted lipid class has
#' more than one chain.
#'
#' Taking PG subclass as an example, intensities of lysoPG fragments
#' (informative for sn1) can be employed to confirm the chains structure
#' (intrules = c("lysopg_sn1/lysopg_sn1")).
#' In this case, the intensity of lysoPG resulting from the loss of the FA chain
#' in sn2 is at least 3 times higher (rates = c("3/1")) than the lysoPG
#' resulting from the loss of the FA chain in sn1.
#'
#' For the intrules argument, "/" will be use to separate the fragments to
#' compare, and "_" will be use to indicate in which list of fragments we
#' need to look for their intensities. This will depend on the chain fragments
#' rules defined previiously.
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
checkIntRules <- function(intrules, rates, intrequired, nchains, combinations,
                          sn1, sn2 = data.frame(), sn3 = data.frame(),
                          sn4 = data.frame()){
  if (nchains == 1){
    passed <- rep(FALSE, length(sn1))
  } else if (length(intrules) == 0 | "Unknown" %in% intrules){
    if (nrow(combinations) > 0){
      passed <- rep(FALSE, nrow(combinations))
      for (c in 1:nrow(combinations)){
        if (nchains == 2){
          if (combinations[c,1] == combinations[c,2]){
            passed[c] <- TRUE
          }
        } else if (nchains == 3){
          if (combinations[c,1] == combinations[c,2] &
              combinations[c,1] == combinations[c,3]){
            passed[c] <- TRUE
          }
        } else if (nchains == 4){
          if (combinations[c,1] == combinations[c,2] &
              combinations[c,1] == combinations[c,3] &
              combinations[c,1] == combinations[c,4]){
            passed[c] <- TRUE
          }
        }
      }
    } else {
      return(F)
    }
  } else if (nchains == 2){
    if (nrow(combinations) == 0){
      return(F)
    } else {
      verified <- matrix(NA,  ncol=length(intrules), nrow=nrow(combinations))
      for (c in 1:nrow(combinations)){
        if (combinations[c,1] == combinations[c,2]){
          verified[c,] <-T
        } else {
          for (i in 1:length(intrules)){
            comp <- unlist(strsplit(intrules[i], "[_/]"))
            list1 <- get(comp[2])
            list2 <- get(comp[4])
            int1 <- as.numeric(list1$int[list1$cb == combinations[c,1] &
                                           list1$db == comp[1]])
            if (length(int1) == 0){int1 <- 0}
            int2 <- as.numeric(list2$int[list2$cb == combinations[c,2] &
                                           list2$db == comp[3]])
            if (length(int2) == 0){int2 <- 0}
            if (length(int1) == 1 & length(int2) == 1){
              if (int1 == 0 & int2 == 0){
                verified[c,i] <- F
              } else {
                if (eval(parse(text=rates[i])) > 1){
                  check <- int1/int2 >= eval(parse(text=rates[i]))
                } else {
                  check <- int1/int2 <= eval(parse(text=rates[i]))
                }
                verified[c,i] <- check
              }
            } else {
              if (eval(parse(text=rates[i])) > 1){
                check <- any(int1[1]/int2 >= eval(parse(text=rates[i])))
              } else {
                check <- any(int1[1]/int2 <= eval(parse(text=rates[i])))
              }
              verified[c,i] <- check
            }
          }
        }
      }
      if (any(intrequired)){ # when there are required fragments, we check those
        passed <- unlist((apply(verified, 1, function(x){
          sum(x)
        }))) >= sum(intrequired)
      } else { # if there isnt any required fragment, any of them will be enough
        passed <- apply(verified, 1, sum) > 0
      }
    }
  } else if (nchains == 3){
    if (nrow(combinations) == 0){
      return(F)
    } else {
      verified <- matrix(NA,  ncol=length(intrules), nrow=nrow(combinations))
      for (c in 1:nrow(combinations)){
        if (combinations[c,1] == combinations[c,2] &
            combinations[c,2] == combinations[c,3]){
          verified[c,] <-T
        } else {
          for (i in 1:length(intrules)){
            comp <- unlist(strsplit(intrules[i], "[_/]"))
            r <- unlist(strsplit(rates[i], "[/]"))
            list1 <- get(comp[2])
            list2 <- get(comp[4])
            list3 <- get(comp[6])
            int1 <- as.numeric(list1$int[list1$cb == combinations[c,1] &
                                           list1$db == comp[1]])
            if (length(int1) == 0){int1 <- 0.1}
            int2 <- as.numeric(list2$int[list2$cb == combinations[c,2] &
                                           list2$db == comp[3]])
            if (length(int2) == 0){int2 <- 0.1}
            int3 <- as.numeric(list3$int[list3$cb == combinations[c,3] &
                                           list3$db == comp[5]])
            if (length(int3) == 0){int3 <- 0.1}
            if (length(int1) == 1 & length(int2) == 1 & length(int3) == 1){
              if (int1 == 0.1 & int2 == 0.1 & int3 == 0.1){
                verified[c,i] <- F
              } else {
                if (as.numeric(r[1])/ as.numeric(r[2]) > 1){
                  check1 <- int1/int2 >= as.numeric(r[1])/ as.numeric(r[2])
                } else {
                  check1 <- int1/int2 <= as.numeric(r[1])/ as.numeric(r[2])
                }
                if (as.numeric(r[2])/ as.numeric(r[3]) > 1){
                  check2 <- int1/int2 >= as.numeric(r[2])/ as.numeric(r[3])
                } else {
                  check2 <- int1/int2 <= as.numeric(r[2])/ as.numeric(r[3])
                }
                check <- (check1 + check2) == 2
                verified[c,i] <- check
              }
            }
          }
        }
      }
      if (any(intrequired)){ # when there are required fragments, we check those
        passed <- unlist((apply(verified, 1, function(x){
          sum(x)
        }))) >= sum(intrequired)
      } else { # if there isnt any required fragment, any of them will be enough
        passed <- apply(verified, 1, sum) > 0
      }
    }
  }
  return(passed)
}

# removeAdducts
#' remove wrong assigned adducts
#'
#' remove wrong assigned adducts.
#'
#' @param df data frame with the input results
#' @param mdiff mass difference expected between adducts
#' @param ppm mass tolerance
#' @param rttol rt tolerance
#'
#'
#' @return Data frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
removeAdducts <- function(df, mdiff, ppm, rttol){
  adducts <- rep(TRUE, nrow(df))
  for (i in 1:nrow(df)){
    adduct <- which(abs((df$m.z - df$m.z[i]) - mdiff)*1000000/df$m.z[i] < ppm)
    if (length(adduct) > 0){
      for (u in 1:length(adduct)){
        if (abs(df[i, "RT"] - df[adduct[u], "RT"]) < rttol &&
            df[adduct[u], "int"] < df[i, "int"]){
          adducts[adduct[u]] <- FALSE
        }
      }
    }
  }
  res <- df[adducts,]
  return(res)
}

# removeDGs
#' remove fragments whose relative intensity to the precursor is less than 10
#' times smaller
#'
#' remove fragments whose relative intensity to the precursor is less than 10
#' times smaller. It is used by idTGpos to reduce false positive annotations.
#'
#' @param sn fragments
#' @param candidates precursors
#'
#'
#' @return Subsetted data frame
#'
#' @keywords internal
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
removeDGs = function(sn, candidates){
  for (c in 1:nrow(candidates)){
    if (nrow(sn[[c]]) > 0){
      sn[[c]] = sn[[c]][as.numeric(sn[[c]]$int)/candidates$int[c] > 0.1,]
    } else {
      sn[[c]] = data.frame()
    }
  }
  return(sn)
}

