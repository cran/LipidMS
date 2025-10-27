## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=T, eval=F-----------------------------------------------------------
# # Install LipidMS
# install.packages("LipidMS", dependencies = c("Depends", "Imports"))
# 
# # load library
# library(LipidMS)
# 

## ----echo=T, eval=F-----------------------------------------------------------
# # load LipidMS library
# library(LipidMS)
# 
# # Example data files can be downloaded from:
# # <https://drive.google.com/drive/folders/1hSYrQBkh-rAA-oiaKqGkrNL7uWQraV75?usp=sharing>
# 
# # get mzXML files from your working directory
# files <- dir()[grepl(".mzXML", dir())]
# 
# # set the processing parameters
# acquisitionmode <- c("DIA", "DDA", "MS", "MS", "MS", "DIA")
# polarity <- "negative"
# dmzagglom <- 15
# drtagglom <- 500
# drtclust <- c(100, 200)
# minpeak <- c(5, 4)
# drtgap <- 5
# drtminpeak <- 15
# drtmaxpeak <- c(100, 200)
# recurs <- c(5, 10)
# sb <- c(3, 2)
# sn <- c(3, 2)
# minint <- c(500, 100)
# weight <- c(2, 3)
# dmzIso <- 5
# drtIso <- 5
# 
# # for single file processing, only samples acquired in DIA or DDA will be processed
# files <- files[acquisitionmode %in% c("DIA", "DDA")]
# acquisitionmode <- acquisitionmode[acquisitionmode %in% c("DIA", "DDA")]
# 
# # run the dataProcessing function to obtain the requires msobjects
# msobjects <- list()
# 
# for (f in 1:length(files)){
#   msobjects[[f]] <- dataProcessing(file = files[f],
#                                    polarity = polarity,
#                                    dmzagglom = dmzagglom,
#                                    drtagglom = drtagglom,
#                                    drtclust = drtclust,
#                                    minpeak = minpeak,
#                                    drtgap = drtgap,
#                                    drtminpeak = drtminpeak,
#                                    drtmaxpeak = drtmaxpeak,
#                                    recurs = recurs,
#                                    sb = sb,
#                                    sn = sn,
#                                    minint = minint,
#                                    weight = weight,
#                                    dmzIso = dmzIso,
#                                    drtIso = drtIso)
# }

## ----echo=T, eval=F-----------------------------------------------------------
# # set annotation parameters
# dmzprecursor <- 5
# dmzproducts <- 10
# rttol <- 5
# coelcutoff <- 0.8
# 
# # Annotate lipids
# ## If polarity is positive
# if (polarity == "positive"){
#   for (m in 1:length(msobjects)){
#     msobjects[[m]] <- idPOS(msobjects[[m]],
#                             ppm_precursor = dmzprecursor,
#                             ppm_products = dmzproducts,
#                             rttol = rttol,
#                             coelCutoff = coelcutoff)
#   }
# }
# 
# ## If polarity is negative
# if (polarity == "negative"){
#   for (m in 1:length(msobjects)){
#     msobjects[[m]] <- idNEG(msobjects[[m]],
#                             ppm_precursor = dmzprecursor,
#                             ppm_products = dmzproducts,
#                             rttol = rttol,
#                             coelCutoff = coelcutoff)
#   }
# }
# 

## ----echo=T, eval=F-----------------------------------------------------------
# # example code for idPEpos function
# pe <- idPEpos(msobject,
#           ppm_precursor = ppm_precursor,
#           ppm_products = ppm_products, rttol = 6,
#           chainfrags_sn1 = c("mg_M+H-H2O", "lysope_M+H-H2O"),
#           chainfrags_sn2 = c("fa_M+H-H2O", "mg_M+H-H2O"),
#           intrules = c("mg_sn1/mg_sn2", "lysope_sn1/lysope_sn2"),
#           rates = c("3/1", "3/1"), intrequired = c(T),
#           dbs = dbs, coelCutoff = 0.8)

## ----echo=T, eval=F-----------------------------------------------------------
# msobject <- idPOS(msobject)
# msobject <- plotLipids(msobject)
# 
# # display the first plot
# msobject$plots[[1]]
# msobject$plots[["yourpeakIDofinterest"]]
# 
# # save all plot to a pdf file
# pdf("plotresults.pdf")
# for (p in 1:length(msobject$plots)){
#     print(msobject$plots[[p]])
#   }
# dev.off()

## ----echo=T, eval=F-----------------------------------------------------------
# # load LipidMS library
# library(LipidMS)
# 
# # Example data files can be downloaded from:
# # <https://drive.google.com/drive/folders/1hSYrQBkh-rAA-oiaKqGkrNL7uWQraV75?usp=sharing>
# 
# # csv file with 3 columns: sample (mzXML file names), acquisitionmode
# # (MS, DIA or DDA) and sampletype (QC, group1, group2, etc.)
# metadata <- read.csv("Matadata.csv", sep=",", dec = ".")
# 
# #==============================================================================#
# # Set processing parameters
# #==============================================================================#
# 
# #===================#
# # Peak-picking
# polarity <- "positive" # 6550 abrir hasta 50 ppm, con el orbi dejar en 15-20 ppm
# dmzagglom <- 15
# drtagglom <- 500
# drtclust <- c(100, 200)
# minpeak <- c(8, 5)
# drtgap <- 5
# drtminpeak <- 15
# drtmaxpeak <- 200
# recurs <- c(5, 10)
# sb <- c(3,2)
# sn <- c(3,2)
# minint <- c(5000, 1000)
# weight <- c(2, 3)
# dmzIso <- 5
# drtIso <- 5
# 
# #==============================================================================#
# # Processing
# #==============================================================================#
# 
# #===================#
# # Peak-picking
# msbatch <- batchdataProcessing(metadata = metadata,
#                                polarity = polarity,
#                                dmzagglom = dmzagglom,
#                                drtagglom = drtagglom,
#                                drtclust = drtclust,
#                                minpeak = minpeak,
#                                drtgap = drtgap,
#                                drtminpeak = drtminpeak,
#                                drtmaxpeak = drtmaxpeak,
#                                recurs = recurs,
#                                sb = sb,
#                                sn = sn,
#                                minint = minint,
#                                weight = weight,
#                                dmzIso = dmzIso,
#                                drtIso = drtIso,
#                                parallel = parallel,
#                                ncores = ncores)
# 
# save(msbatch, file="msbatch.rda.gz", compress = TRUE)

## ----echo=T, eval=F-----------------------------------------------------------
# #==============================================================================#
# # Set parameters
# #==============================================================================#
# 
# #===================#
# # Batch processing
# dmzalign <- 5
# drtalign <- 30
# span <- 0.4
# minsamplesfracalign <- 0.75
# dmzgroup <- 5
# drtagglomgroup <- 30
# drtgroup <- 15
# minsamplesfracgroup <- 0.25
# parallel <- TRUE
# ncores <- 2
# 
# #==============================================================================#
# # Processing
# #==============================================================================#
# 
# #===================#
# # Alignment
# msbatch <- alignmsbatch(msbatch, dmz = dmzalign, drt = drtalign, span = span,
#                         minsamplesfrac = minsamplesfracalign,
#                         parallel = parallel, ncores = ncores)
# 
# # rt deviation plot
# rtdevplot(msbatch)
# rtdevplot(msbatch, colorbygroup = FALSE)
# 
# # tic plot
# plotticmsbatch(msbatch)
# plotticmsbatch(msbatch, colorbygroup = FALSE)
# 
# #===================#
# # Grouping
# msbatch <- groupmsbatch(msbatch, dmz = dmzgroup, drtagglom = drtagglomgroup,
#                         drt = drtgroup, minsamplesfrac = minsamplesfracgroup,
#                         parallel = parallel, ncores = ncores)
# 
# 
# #===================#
# # Fill missing peaks
# msbatch <- fillpeaksmsbatch(msbatch)
# 
# 
# # Now we have a data matrix with all samples and features.
# View(msbatch$features)

## ----echo=T, eval=F-----------------------------------------------------------
# #==============================================================================#
# # Lipid annotation
# #==============================================================================#
# 
# #===================#
# # Lipid Annotation
# msbatch <- annotatemsbatch(msbatch)
# 
# 
# # Make plots for identified lipids
# for (m in 1:length(msbatch$msobjects)){
#   if (msbatch$msobjects[[m]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
#     msbatch$msobjects[[m]] <- plotLipids(msbatch$msobjects[[m]])
#   }
# }
# 
# print(msbatch$msobjects[[1]]$annotation$plots[[1]])
# 

## ----echo=T, eval=F-----------------------------------------------------------
# #===================#
# # features
# peaklist <- msbatch$features
# peaklistNoIso <- peaklist[peaklist$isotope %in% c("", "[M+0]"),]
# 
# View(peaklistNoIso)
# 
# write.csv(peaklist, file="peaklist.csv")
# write.csv(peaklistNoIso, file="peaklistNoIso.csv")
# 
# #===================#
# # annotations
# 
# # results
# for (i in 1:length(msbatch$msobjects)){
#   if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
#     fileName <- gsub(".mzXML", "_summaryResults.csv" ,
#                      msbatch$msobjects[[i]]$metaData$generalMetadata$file[i])
#     write.csv(msbatch$msobjects[[i]]$annotation$results, fileName, row.names = FALSE)
#   }
# }
# 
# # Annotated Peaklists
# for (i in 1:length(msbatch$msobjects)){
#   if (msbatch$msobjects[[i]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
#     fileName <- gsub(".mzXML", "_annotatedPeaklist.csv" ,
#                      msbatch$msobjects[[i]]$metaData$generalMetadata$file[i])
#     write.csv(msbatch$msobjects[[i]]$annotation$annotatedPeaklist, fileName, row.names = FALSE)
#   }
# }
# 
# #===================#
# # plots
# 
# pdf("RTdevplot.pdf")
# rtdevplot(msbatch)
# rtdevplot(msbatch, colorbygroup = FALSE)
# dev.off()
# 
# pdf("TIC.pdf", height = 7, width = 10)
# plotticmsbatch(msbatch)
# plotticmsbatch(msbatch, colorbygroup = FALSE)
# dev.off()
# 
# # lipid id plots
# for (s in 1:length(msbatch$msobjects)){
#   if (msbatch$msobjects[[s]]$metaData$generalMetadata$acquisitionmode %in% c("DIA", "DDA")){
#     print(s)
#     if (msbatch$msobjects[[s]]$metaData$generalMetadata$acquisitionmode == "DIA"){
#       height <- 7
#     } else {
#       height <- 9
#     }
#     pdf(file = gsub(".mzXML", "_plots.pdf", msbatch$msobjects[[s]]$metaData$generalMetadata$file),
#         width = 8, height = height)
#     for ( pl in 1:length(msbatch$msobjects[[s]]$annotation$plots)){
#       print(msbatch$msobjects[[s]]$annotation$plots[[pl]])
#     }
#     dev.off()
#   }
# }

## ----echo=T, eval=F-----------------------------------------------------------
# # To run LipidMS shiny app execute:
# LipidMSapp()

## ----echo=T, eval=F-----------------------------------------------------------
# fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
# "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "20:0", "20:1", "20:2",
# "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4", "22:5",
# "22:6", "24:0", "24:1", "26:0")
# sph <- c("16:0", "16:1", "18:0", "18:1")
# dbs <- createLipidDB(lipid = "all", chains = fas, chains2 = sph)
# 
# # to use for identification function two additional data frames need to be added
# dbs$adductsTable <- LipidMS::adductsTable
# dbs$nlsphdb <- LipidMS::nlsphdb

## ----echo=T, eval=F-----------------------------------------------------------
# fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#          "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "19:0", "20:0", "20:1",
#          "20:2", "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4",
#          "22:5", "22:6", "24:0", "24:1", "26:0")
# newfadb <- createLipidDB(lipid = "FA", chains = fas)
# dbs <- assignDB() # This function loads all DBs required
# dbs$fadb <- newfadb$fadb # Then, you can modify some of these DBs

## ----echo=T, eval=F-----------------------------------------------------------
# adductsTable <- LipidMS::adductsTable
# adductsTable <- data.frame(adduct = c(adductsTable$adduct, "M+X"),
#                           mdiff = c(adductsTable$mdiff, 52.65),
#                           charge = c(adductsTable$charge, 1),
#                           n = c(adductsTable$n, 1),
#                           stringsAsFactors = F)

## ----echo=T, eval=F-----------------------------------------------------------
# # The new adductsTable has to be also uploaded in the dbs list.
# dbs <- assignDB()
# dbs$adductsTable <- adductsTable
# 
# idPCpos(msobject = LipidMSdata2::msobjectDIApos,
#         adducts = c("M+H", "M+Na", "M+X"), dbs = dbs)

