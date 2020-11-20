# createLipidDB
#' Customizable lipid DBs creator
#'
#' It allows to create easy-customizable lipid DBs for annotation with LipidMS
#' package.
#'
#' @param lipid character value indicating the class of lipid. See Details.
#' @param chains character vector indicating the FA chains to be employed
#' @param chains2 character vector containing the sphingoid bases to be employed
#' if required.
#'
#' @return List with the requested dbs (data frames)
#'
#' @details \code{lipidClass} argument needs to be one of the following
#' character values: "Cer", "CerP", "GlcCer", "SM", "Carnitine", "CE", "FA",
#' "HFA", "Sph" (sphingoid bases), "SphP", "MG", "LPA", , "LPC",
#' "LPE", "LPG", "LPI", "LPS", "FAHFA", "DG", "PC", "PE", "PG", "PI", "PS",
#' "PA", "TG", "CL" or "all".
#'
#' @examples
#' fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#' "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "20:0", "20:1", "20:2",
#' "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4", "22:5",
#' "22:6", "24:0", "24:1", "26:0")
#' sph <- c("16:0", "16:1", "18:0", "18:1")
#' newdb <- createLipidDB(lipid = "PC", chains = fas, chains2 = sph)
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
createLipidDB <- function(lipid, chains, chains2){
  customizedDataSets <- list()
  if (sum(lipid == "Cer" || lipid == "CerP" || lipid == "GlcCer" ||
          lipid == "SM") == 1){
    db <- dbSphingolipids(chains = chains, chains2 = chains2, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    if (lipid == "Cer"){
      customizedDataSets[["cerdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "CerP"){
      customizedDataSets[["cerPdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "GlcCer"){
      customizedDataSets[["glccerdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "SM"){
      customizedDataSets[["smdb"]] <- data.frame(formula=db$formula,
                                                 total=db$total, Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
  } else if (sum(lipid == "FA" || lipid == "HFA" || lipid == "Carnitine" ||
                 lipid == "LPA" || lipid == "LPE" || lipid == "LPG" ||
                 lipid == "LPI" || lipid == "LPS" || lipid == "LPC" ||
                 lipid == "MG" || lipid == "CE" || lipid == "Sph" ||
                 lipid == "SphP" || lipid == "LPEo" || lipid == "LPAo" ||
                 lipid == "LPCp" || lipid == "LPCo" || lipid == "LPEp") == 1){
    if (lipid %in% c("Sph", "SphP")){
      db <- dbOneChain(chains = chains2, lipid = lipid)
      db <- data.frame(formula=db$formula, total=db$total,
                       Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                       stringsAsFactors = F)
    } else {
      db <- dbOneChain(chains = chains, lipid = lipid)
      db <- data.frame(formula=db$formula, total=db$total,
                       Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                       stringsAsFactors = F)
    }
    if (lipid == "FA"){
      customizedDataSets[["fadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "HFA"){
      customizedDataSets[["hfadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "Carnitine"){
      customizedDataSets[["carnitinedb"]] <- data.frame(formula=db$formula, total=db$total,
                                                        Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPA"){
      customizedDataSets[["lysopadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPE"){
      customizedDataSets[["lysopedb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPG"){
      customizedDataSets[["lysopgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPI"){
      customizedDataSets[["lysopidb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPS"){
      customizedDataSets[["lysopsdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPC"){
      customizedDataSets[["lysopcdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                     Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "MG"){
      customizedDataSets[["mgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "CE"){
      customizedDataSets[["CEdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "Sph"){
      customizedDataSets[["sphdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "SphP"){
      customizedDataSets[["sphPdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPEo"){
      customizedDataSets[["lysopeodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPAo"){
      customizedDataSets[["lysopaodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPCp"){
      customizedDataSets[["lysopcpdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPCo"){
      customizedDataSets[["lysopcodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "LPEp"){
      customizedDataSets[["lysopepdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
  } else if (sum(lipid == "DG" || lipid == "PC" || lipid == "PE" ||
                 lipid == "PG" || lipid == "PI" || lipid == "PS" || lipid == "PIP" ||
                 lipid == "PIP2" ||lipid == "PIP3" || lipid == "FAHFA" ||
                 lipid == "PA") == 1 || lipid == "PEo" || lipid == "PCp" ||
             lipid == "PCo" || lipid == "PEp"){
    db <- dbTwoChains(chains = chains, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    if (lipid == "FAHFA"){
      customizedDataSets[["fahfadb"]] <- data.frame(formula=db$formula, total=db$total,
                                                    Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "DG"){
      customizedDataSets[["dgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PE"){
      customizedDataSets[["pedb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PG"){
      customizedDataSets[["pgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PI"){
      customizedDataSets[["pidb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PIP"){
      customizedDataSets[["pipdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                  Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PIP2"){
      customizedDataSets[["pip2db"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PIP3"){
      customizedDataSets[["pip3db"]] <- data.frame(formula=db$formula, total=db$total,
                                                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PS"){
      customizedDataSets[["psdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PC"){
      customizedDataSets[["pcdb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PA"){
      customizedDataSets[["padb"]] <- data.frame(formula=db$formula, total=db$total,
                                                 Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PEo"){
      customizedDataSets[["peodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PCp"){
      customizedDataSets[["pcpdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PCo"){
      customizedDataSets[["pcodb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
    if (lipid == "PEp"){
      customizedDataSets[["pepdb"]] <-
        data.frame(formula=db$formula, total=db$total,
                   Mass=as.numeric(db$Mass), stringsAsFactors = F)
    }
  } else if (lipid == "TG") {
    db <- dbThreeChains(chains = chains, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    customizedDataSets[["tgdb"]] <- data.frame(formula=db$formula, total=db$total,
                                               Mass=as.numeric(db$Mass), stringsAsFactors = F)
  } else if (lipid == "CL") {
    db <- dbFourChains(chains = chains, lipid = lipid)
    db <- data.frame(formula=db$formula, total=db$total,
                     Mass=as.numeric(db$Mass), ID = paste(lipid, "(", db$total, ")", sep=""),
                     stringsAsFactors = F)
    customizedDataSets[["cldb"]] <- data.frame(formula=db$formula, total=db$total,
                                               Mass=as.numeric(db$Mass), stringsAsFactors = F)
  } else if (lipid == "all"){
    ceramides <- dbSphingolipids(chains = chains, chains2 = chains2,
                                 lipid = "Cer")
    ceramides <- data.frame(formula=ceramides$formula, total=ceramides$total,
                            Mass=as.numeric(ceramides$Mass), ID = paste("Cer(", ceramides$total,
                                                                        ")", sep=""), stringsAsFactors = F)
    ceramidesP <- dbSphingolipids(chains = chains, chains2 = chains2,
                                  lipid = "CerP")
    ceramidesP <- data.frame(formula=ceramidesP$formula, total=ceramidesP$total,
                             Mass=as.numeric(ceramidesP$Mass), ID = paste("CerP(", ceramidesP$total,
                                                                          ")", sep=""), stringsAsFactors = F)
    glccer <- dbSphingolipids(chains = chains, chains2 = chains2,
                              lipid = "GlcCer")
    glccer <- data.frame(formula=glccer$formula, total=glccer$total,
                         Mass=as.numeric(glccer$Mass), ID = paste("GlcCer(", glccer$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    sm <- dbSphingolipids(chains = chains, chains2 = chains2, lipid = "SM")
    sm <- data.frame(formula=sm$formula, total=sm$total,
                     Mass=as.numeric(sm$Mass), ID = paste("SM(", sm$total,
                                                          ")", sep=""), stringsAsFactors = F)
    fa <- dbOneChain(chains = chains, lipid = "FA")
    fa <- data.frame(formula=fa$formula, total=fa$total,
                     Mass=as.numeric(fa$Mass), ID = paste("FA(", fa$total,
                                                          ")", sep=""), stringsAsFactors = F)
    hfa <- dbOneChain(chains = chains, lipid = "HFA")
    hfa <- data.frame(formula=hfa$formula, total=hfa$total,
                      Mass=as.numeric(hfa$Mass), ID = paste("HFA(", hfa$total,
                                                            ")", sep=""), stringsAsFactors = F)
    carnitine <- dbOneChain(chains = chains, lipid = "Carnitine")
    carnitine <- data.frame(formula=carnitine$formula, total=carnitine$total,
                            Mass=as.numeric(carnitine$Mass), ID = paste("Carnitine(", carnitine$total,
                                                                        ")", sep=""), stringsAsFactors = F)
    CE <- dbOneChain(chains = chains, lipid = "CE")
    CE <- data.frame(formula=CE$formula, total=CE$total,
                     Mass=as.numeric(CE$Mass), ID = paste("CE(", CE$total,
                                                          ")", sep=""), stringsAsFactors = F)
    mg <- dbOneChain(chains = chains, lipid = "MG")
    mg <- data.frame(formula=mg$formula, total=mg$total,
                     Mass=as.numeric(mg$Mass), ID = paste("MG(", mg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    sph <- dbOneChain(chains = chains2, lipid = "Sph")
    sph <- data.frame(formula=sph$formula, total=sph$total,
                      Mass=as.numeric(sph$Mass), ID = paste("Sph(", sph$total,
                                                            ")", sep=""), stringsAsFactors = F)
    sphP <- dbOneChain(chains = chains2, lipid = "SphP")
    sphP <- data.frame(formula=sphP$formula, total=sphP$total,
                       Mass=as.numeric(sphP$Mass), ID = paste("SphP(", sphP$total,
                                                              ")", sep=""), stringsAsFactors = F)
    lysopc <- dbOneChain(chains = chains, lipid = "LPC")
    lysopc <- data.frame(formula=lysopc$formula, total=lysopc$total,
                         Mass=as.numeric(lysopc$Mass), ID = paste("LPC(", lysopc$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysope <- dbOneChain(chains = chains, lipid = "LPE")
    lysope <- data.frame(formula=lysope$formula, total=lysope$total,
                         Mass=as.numeric(lysope$Mass), ID = paste("LPE(", lysope$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysopg <- dbOneChain(chains = chains, lipid = "LPG")
    lysopg <- data.frame(formula=lysopg$formula, total=lysopg$total,
                         Mass=as.numeric(lysopg$Mass), ID = paste("LPG(", lysopg$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysopi <- dbOneChain(chains = chains, lipid = "LPI")
    lysopi <- data.frame(formula=lysopi$formula, total=lysopi$total,
                         Mass=as.numeric(lysopi$Mass), ID = paste("LPI(", lysopi$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysops <- dbOneChain(chains = chains, lipid = "LPS")
    lysops <- data.frame(formula=lysops$formula, total=lysops$total,
                         Mass=as.numeric(lysops$Mass), ID = paste("LPS(", lysops$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    lysopa <- dbOneChain(chains = chains, lipid = "LPA")
    lysopa <- data.frame(formula=lysopa$formula, total=lysopa$total,
                         Mass=as.numeric(lysopa$Mass), ID = paste("LPA(", lysopa$total,
                                                                  ")", sep=""), stringsAsFactors = F)
    pc <- dbTwoChains(chains = chains, lipid = "PC")
    pc <- data.frame(formula=pc$formula, total=pc$total,
                     Mass=as.numeric(pc$Mass), ID = paste("PC(", pc$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pe <- dbTwoChains(chains = chains, lipid = "PE")
    pe <- data.frame(formula=pe$formula, total=pe$total,
                     Mass=as.numeric(pe$Mass), ID = paste("PE(", pe$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pg <- dbTwoChains(chains = chains, lipid = "PG")
    pg <- data.frame(formula=pg$formula, total=pg$total,
                     Mass=as.numeric(pg$Mass), ID = paste("PG(", pg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pi <- dbTwoChains(chains = chains, lipid = "PI")
    pi <- data.frame(formula=pi$formula, total=pi$total,
                     Mass=as.numeric(pi$Mass), ID = paste("PI(", pi$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pip <- dbTwoChains(chains = chains, lipid = "PIP")
    pip <- data.frame(formula=pip$formula, total=pip$total,
                      Mass=as.numeric(pip$Mass), ID = paste("PIP(", pip$total,
                                                            ")", sep=""), stringsAsFactors = F)
    pip2 <- dbTwoChains(chains = chains, lipid = "PIP2")
    pip2 <- data.frame(formula=pip2$formula, total=pip2$total,
                       Mass=as.numeric(pip2$Mass), ID = paste("PIP2(", pip2$total,
                                                              ")", sep=""), stringsAsFactors = F)
    pip3 <- dbTwoChains(chains = chains, lipid = "PIP3")
    pip3 <- data.frame(formula=pip3$formula, total=pip3$total,
                       Mass=as.numeric(pip3$Mass), ID = paste("PIP3(", pip3$total,
                                                              ")", sep=""), stringsAsFactors = F)
    ps <- dbTwoChains(chains = chains, lipid = "PS")
    ps <- data.frame(formula=ps$formula, total=ps$total,
                     Mass=as.numeric(ps$Mass), ID = paste("PS(", ps$total,
                                                          ")", sep=""), stringsAsFactors = F)
    pa <- dbTwoChains(chains = chains, lipid = "PA")
    pa <- data.frame(formula=pa$formula, total=pa$total,
                     Mass=as.numeric(pa$Mass), ID = paste("PA(", pa$total,
                                                          ")", sep=""), stringsAsFactors = F)
    fahfa <- dbTwoChains(chains = chains, lipid = "FAHFA")
    fahfa <- data.frame(formula=fahfa$formula, total=fahfa$total,
                        Mass=as.numeric(fahfa$Mass), ID = paste("FAHFA(", fahfa$total,
                                                                ")", sep=""), stringsAsFactors = F)
    dg <- dbTwoChains(chains = chains, lipid = "DG")
    dg <- data.frame(formula=dg$formula, total=dg$total,
                     Mass=as.numeric(dg$Mass), ID = paste("DG(", dg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    tg <- dbThreeChains(chains = chains, lipid = "TG")
    tg <- data.frame(formula=tg$formula, total=tg$total,
                     Mass=as.numeric(tg$Mass), ID = paste("TG(", tg$total,
                                                          ")", sep=""), stringsAsFactors = F)
    cl <- dbFourChains(chains = chains, lipid = "CL")
    cl <- data.frame(formula=cl$formula, total=cl$total,
                     Mass=as.numeric(cl$Mass), ID = paste("CL(", cl$total,
                                                          ")", sep=""), stringsAsFactors = F)
    peo <- dbTwoChains(chains = chains, lipid = "PEo")
    peo <- data.frame(formula=peo$formula, total=peo$total,
                      Mass=as.numeric(peo$Mass),
                      ID = paste("PEo(", peo$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopeo <- dbOneChain(chains = chains, lipid = "LPEo")
    lysopeo <- data.frame(formula=lysopeo$formula, total=lysopeo$total,
                          Mass=as.numeric(lysopeo$Mass),
                          ID = paste("LPEo(", pe$total, ")", sep=""),
                          stringsAsFactors = F)
    lysopao <- dbOneChain(chains = chains, lipid = "LPAo")
    lysopao <- data.frame(formula=lysopao$formula, total=lysopao$total,
                          Mass=as.numeric(lysopao$Mass),
                          ID = paste("LPAo(", pe$total, ")", sep=""),
                          stringsAsFactors = F)
    pcp <- dbTwoChains(chains = chains, lipid = "PCp")
    pcp <- data.frame(formula=pcp$formula, total=pcp$total,
                      Mass=as.numeric(pcp$Mass),
                      ID = paste("PCp(", pcp$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopcp <- dbOneChain(chains = chains, lipid = "LPCp")
    lysopcp <- data.frame(formula=lysopcp$formula, total=lysopcp$total,
                          Mass=as.numeric(lysopcp$Mass),
                          ID = paste("LPCp(", lysopcp$total, ")", sep=""),
                          stringsAsFactors = F)
    pco <- dbTwoChains(chains = chains, lipid = "PCo")
    pco <- data.frame(formula=pco$formula, total=pco$total,
                      Mass=as.numeric(pco$Mass),
                      ID = paste("PCo(", pco$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopco <- dbOneChain(chains = chains, lipid = "LPCo")
    lysopco <- data.frame(formula=lysopco$formula, total=lysopco$total,
                          Mass=as.numeric(lysopco$Mass),
                          ID = paste("LPCo(", lysopco$total, ")", sep=""),
                          stringsAsFactors = F)
    pep <- dbTwoChains(chains = chains, lipid = "PEp")
    pep <- data.frame(formula=pep$formula, total=pep$total,
                      Mass=as.numeric(pep$Mass),
                      ID = paste("PEp(", pco$total, ")", sep=""),
                      stringsAsFactors = F)
    lysopep <- dbOneChain(chains = chains, lipid = "LPEp")
    lysopep <- data.frame(formula=lysopep$formula, total=lysopep$total,
                          Mass=as.numeric(lysopep$Mass),
                          ID = paste("LPEp(", lysopep$total, ")", sep=""),
                          stringsAsFactors = F)
    customizedDataSets[["cerdb"]] <- ceramides
    customizedDataSets[["cerPdb"]] <- ceramidesP
    customizedDataSets[["glccerdb"]] <- glccer
    customizedDataSets[["smdb"]] <- sm
    customizedDataSets[["fadb"]] <- fa
    customizedDataSets[["hfadb"]] <- hfa
    customizedDataSets[["carnitinedb"]] <- carnitine
    customizedDataSets[["lysopadb"]] <- lysopa
    customizedDataSets[["lysopedb"]] <- lysope
    customizedDataSets[["lysopgdb"]] <- lysopg
    customizedDataSets[["lysopidb"]] <- lysopi
    customizedDataSets[["lysopsdb"]] <- lysops
    customizedDataSets[["lysopcdb"]] <- lysopc
    customizedDataSets[["mgdb"]] <- mg
    customizedDataSets[["CEdb"]] <- CE
    customizedDataSets[["sphdb"]] <- sph
    customizedDataSets[["sphPdb"]] <- sphP
    customizedDataSets[["fahfadb"]] <- fahfa
    customizedDataSets[["pedb"]] <- pe
    customizedDataSets[["pgdb"]] <- pg
    customizedDataSets[["pidb"]] <- pi
    customizedDataSets[["pipdb"]] <- pip
    customizedDataSets[["pip2db"]] <- pip2
    customizedDataSets[["pip3db"]] <- pip3
    customizedDataSets[["psdb"]] <- ps
    customizedDataSets[["padb"]] <- pa
    customizedDataSets[["pcdb"]] <- pc
    customizedDataSets[["dgdb"]] <- dg
    customizedDataSets[["tgdb"]] <- tg
    customizedDataSets[["cldb"]] <- cl
    customizedDataSets[["peodb"]] <- peo
    customizedDataSets[["lysopeodb"]] <- lysopeo
    customizedDataSets[["lysopaodb"]] <- lysopao
    customizedDataSets[["pcpdb"]] <- pcp
    customizedDataSets[["lysopcpdb"]] <- lysopcp
    customizedDataSets[["pcodb"]] <- pco
    customizedDataSets[["lysopcodb"]] <- lysopco
    customizedDataSets[["pepdb"]] <- pep
    customizedDataSets[["lysopepdb"]] <- lysopep
  }
  return(customizedDataSets)
}
