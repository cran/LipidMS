## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=T, eval=F----------------------------------------------------------
#  # Install LipidMS and LipidMSdata2
#  install.packages("LipidMS", dependencies = c("Depends", "Imports"))
#  devtools::install_github("maialba3/LipidMSdata2")
#  
#  # load LipidMS library
#  library(LipidMS)
#  
#  # load example datasets
#  library(LipidMSdata2)

## ---- echo=T, eval=F----------------------------------------------------------
#  # load LipidMS library
#  library(LipidMS)
#  
#  # get mzXML files from your working directory
#  files <- dir()[grepl(".mzXML", dir())]
#  
#  # set the processing parameters
#  acquisitionmode <- "DIA"
#  polarity <- "negative"
#  dmzagglom <- 5
#  drtagglom <- 25
#  drtclust <- 25
#  minpeak <- c(4, 3)
#  drtgap <- 5
#  drtminpeak <- 20
#  drtmaxpeak <- 200
#  recurs <- 5
#  sb <- 2
#  sn <- 2
#  minint <- c(500, 100)
#  weight <- c(2, 3)
#  dmzIso <- 5
#  drtIso <- 5
#  
#  # run the dataProcessing function to obtain the requires msobjects
#  msobjects <- list()
#  
#  for (f in 1:length(files)){
#    msobjects[[f]] <- dataProcessing(file = files[f],
#                                     polarity = polarity,
#                                     dmzagglom = dmzagglom,
#                                     drtagglom = drtagglom,
#                                     drtclust = drtclust,
#                                     minpeak = minpeak,
#                                     drtgap = drtgap,
#                                     drtminpeak = drtminpeak,
#                                     drtmaxpeak = drtmaxpeak,
#                                     recurs = recurs,
#                                     sb = sb,
#                                     sn = sn,
#                                     minint = minint,
#                                     weight = weight,
#                                     dmzIso = dmzIso,
#                                     drtIso = drtIso)
#  }

## ---- echo=T, eval=F----------------------------------------------------------
#  # Use the example datasets (msobjects)
#  devtools::install_github("maialba3/LipidMSdata2")
#  library(LipidMSdata2)
#  
#  # set annotation parameters
#  dmzprecursor <- 5
#  dmzproducts <- 10
#  rttol <- 5
#  coelcutoff <- 0.8
#  
#  # Annotate lipids
#  annotated_msobject <- idPOS(LipidMSdata2::msobjectDIApos,
#                              ppm_precursor = dmzprecursor,
#                              ppm_products = dmzproducts,
#                              rttol = rttol,
#                              coelCutoff = coelcutoff)
#  
#  
#  # Try with your own data
#  ## If polarity is positive
#  if (polarity == "positive"){
#    for (m in 1:length(msobjects)){
#      msobjects[[m]] <- idPOS(msobjects[[m]],
#                              ppm_precursor = dmzprecursor,
#                              ppm_products = dmzproducts,
#                              rttol = rttol,
#                              coelCutoff = coelcutoff)
#    }
#  }
#  
#  ## If polarity is negative
#  if (polarity == "negative"){
#    for (m in 1:length(msobjects)){
#      msobjects[[m]] <- idNEG(msobjects[[m]],
#                              ppm_precursor = dmzprecursor,
#                              ppm_products = dmzproducts,
#                              rttol = rttol,
#                              coelCutoff = coelcutoff)
#    }
#  }
#  

## ---- echo=T, eval=F----------------------------------------------------------
#  # example code for idPEpos function
#  pe <- idPEpos(msobject = LipidMSdata2::msobjectDIApos,
#            ppm_precursor = ppm_precursor,
#            ppm_products = ppm_products, rttol = 6,
#            chainfrags_sn1 = c("mg_M+H-H2O", "lysope_M+H-H2O"),
#            chainfrags_sn2 = c("fa_M+H-H2O", "mg_M+H-H2O"),
#            intrules = c("mg_sn1/mg_sn2", "lysope_sn1/lysope_sn2"),
#            rates = c("3/1", "3/1"), intrequired = c(T),
#            dbs = dbs, coelCutoff = 0.8)
#  
#  # additional information about how to change rules is given in the documentation
#  # of the following functions: chainFrags , checkClass, checkIntensityRules,
#  # coelutingFrags, ddaFrags, combineChains and organizeResults. These functions
#  # could be also empoyed to build customized identification functions.

## ---- echo=T, eval=F----------------------------------------------------------
#  msobject <- idPOS(LipidMSdata2::msobjectDIApos)
#  msobject <- plotLipids(msobject)
#  
#  # display the first plot
#  msobject$plots[[1]]
#  msobject$plots[["yourpeakIDofinterest"]]
#  
#  # save all plot to a pdf file
#  pdf("plotresults.pdf")
#  for (p in 1:length(msobject$plots)){
#      print(msobject$plots[[p]])
#    }
#  dev.off()

## ---- echo=T, eval=F----------------------------------------------------------
#  # To run LipidMS shiny app execute:
#  LipidMSapp()

## ---- echo=T, eval=F----------------------------------------------------------
#  fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#  "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "20:0", "20:1", "20:2",
#  "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4", "22:5",
#  "22:6", "24:0", "24:1", "26:0")
#  sph <- c("16:0", "16:1", "18:0", "18:1")
#  dbs <- createLipidDB(lipid = "all", chains = fas, chains2 = sph)
#  
#  # to use for identification function two additional data frames need to be added
#  dbs$adductsTable <- LipidMS::adductsTable
#  dbs$nlsphdb <- LipidMS::nlsphdb

## ---- echo=T, eval=F----------------------------------------------------------
#  fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#           "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "19:0", "20:0", "20:1",
#           "20:2", "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4",
#           "22:5", "22:6", "24:0", "24:1", "26:0")
#  newfadb <- createLipidDB(lipid = "FA", chains = fas)
#  dbs <- assignDB() # This function loads all DBs required
#  dbs$fadb <- newfadb$fadb # Then, you can modify some of these DBs

## ---- echo=T, eval=F----------------------------------------------------------
#  adductsTable <- LipidMS::adductsTable
#  adductsTable <- data.frame(adduct = c(adductsTable$adduct, "M+X"),
#                            mdiff = c(adductsTable$mdiff, 52.65),
#                            charge = c(adductsTable$charge, 1),
#                            n = c(adductsTable$n, 1),
#                            stringsAsFactors = F)

## ---- echo=T, eval=F----------------------------------------------------------
#  # The new adductsTable has to be also uploaded in the dbs list.
#  dbs <- assignDB()
#  dbs$adductsTable <- adductsTable
#  
#  idPCpos(msobject = LipidMSdata2::msobjectDIApos,
#          adducts = c("M+H", "M+Na", "M+X"), dbs = dbs)

