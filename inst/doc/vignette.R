## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.cap = TRUE
)

## ---- echo=T, eval=F-----------------------------------------------------
#  library(LipidMS)
#  sepByCE(input = "mix_pos.mzXML",  output = "mix_pos_sep")
#  
#  # to convert a batch of files you can use the following code:
#  files <- dir()[grepl(".mzXML", dir())]
#  outputs <- unlist(lapply(sapply(files, strsplit, ".mzXML"), "[[", 1))
#  outputs <- paste(outputs, "_sep", sep="")
#  mapply(sepByCE, files, outputs)

## ---- echo=T, eval=F-----------------------------------------------------
#  ## 1) remove lockspray function
#  
#  # first, you will need to set your working directory where all you raw files are saved. Then run the following code:
#  
#  folders <- dir()[grepl("raw", dir())]
#  
#  for (f in folders){
#    files <- dir(f)
#    unlink(paste(f, files[grep("FUNC003", files)], sep = "/"))
#  }
#  
#  
#  ## 2) convert to mzXML format using MSConvert
#  
#  ## 3) separate mzXML files
#  
#  files <- dir()
#  files <- files[grep(".mzXML", files)]
#  outputs <- paste(unlist(sapply(files, strsplit, ".mzXML")), "_sep", sep="")
#  
#  for (f in 1:length(files)){
#    nCE <- 2  ## change if you have more than 2 remaining functions after removing lockspray
#    lines <- readLines(files[f])
#    scans <- grep("    <scan num", lines)
#    length(scans)
#    runs <- grep("</msRun>", lines)
#    for (i in 1:nCE){
#      pos <- seq(i,length(scans), nCE)
#      lines2write <- c(1:(scans[1]-1))
#      for (x in pos[1:length(pos)-1]){
#        if (x == scans[length(scans)]){
#          lines2write <- append(lines2write, c(scans[x]:c(runs-1)))
#        } else {
#          lines2write <- append(lines2write, c(scans[x]:(scans[x+1]-1)))
#        }
#      }
#      lines2write <- append(lines2write, c(runs:length(lines)))
#      write(lines[lines2write], file=paste(c(outputs[f], as.character(i), ".mzXML"),
#                                           collapse=""))
#    }
#  }
#  

## ---- echo=T, eval=F-----------------------------------------------------
#  ms1 <- dir()[grepl("fullMS.mzXML", dir())] #fullMS or any other nomenclature employed for MS1
#  ms2 <- dir()[grepl("sep40.mzXML"), dir()] #sep40 or any other nomenclature employed for MS2
#  data_ms1 <- sapply(ms1, dataProcessing, msLevel = 1, polarity = "positive")
#  data_ms2 <- sapply(ms2, dataProcessing, msLevel = 2, polarity = "positive")
#  head(data_ms1[[1]]$peaklist)
#  head(data_ms1[[1]]$rawScans)

## ---- echo=T, eval=F-----------------------------------------------------
#  pos_res <- idPOS(MS1 = data_ms1[[1]], MSMS1 = data_ms2[[1]], ppm_precursor = 10,
#                   ppm_products = 10, rttol = 10, coelCutoff = 0.8)

## ---- echo=T, eval=F-----------------------------------------------------
#  MS1 <- data_ms1[[1]] # CE = 0
#  MSMS1 <- data_ms2[[1]] # CE > 0
#  ppm_precursor <- 10
#  ppm_products <- 10
#  rttol <- 10
#  dbs <- assignDB()
#  
#  # example code for idPEpos function
#  pe <- idPEpos(MS1 = MS1, MSMS1 = MSMS1,
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
#  # coelutingFrag, combineChains and organizeResults. These functions could be also
#  # empoyed to build customized identification functions.

## ---- echo=T, eval=F-----------------------------------------------------
#  MS1 <- data_ms1[[1]] # CE = 0
#  MSMS1 <- data_ms2[[1]] # CE > 0
#  ppm_precursor <- 10
#  ppm_products <- 10
#  rttol <- 10
#  
#  # example to customize several id functions
#  results <- vector()
#  results <- rbind(results, idLPCpos(MS1 = MS1, MSMS1 = MSMS1,
#            ppm_precursor = ppm_precursor,
#            ppm_products = ppm_products, rttol = rttol)$results)
#  results <- rbind(results, idLPEpos(MS1 = MS1, MSMS1 = MSMS1,
#            ppm_precursor = ppm_precursor,
#            ppm_products = ppm_products, rttol = rttol)$results)
#  # in this case you should add the rest of functions for annotation in ESI+
#  
#  # Once you have the complete results table, you can cross it with the MS1 original peak table.
#  annotatedPeaklist <- crossTables(MS1, results, ppm_precursor, rttol)
#  
#  
#  # To see results:
#  View(results)
#  View(annotatedPeaklist)
#  
#  # if you want to see just the features annotated by LipidMS:
#  View(annotatedPeaklist[annotatedPeaklist$LipidMS_id != "",])

## ---- echo=T, eval=F-----------------------------------------------------
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

## ---- echo=T, eval=F-----------------------------------------------------
#  fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
#           "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "19:0", "20:0", "20:1",
#           "20:2", "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4",
#           "22:5", "22:6", "24:0", "24:1", "26:0")
#  newfadb <- createLipidDB(lipid = "FA", chains = fas)
#  dbs <- assignDB() # This function loads all DBs required
#  dbs$fadb <- newfadb$fadb # Then, you can modify some of these DBs

## ---- echo=T, eval=F-----------------------------------------------------
#  adductsTable <- LipidMS::adductsTable
#  adductsTable <- data.frame(adduct = c(adductsTable$adduct, "M+X"),
#                            mdiff = c(adductsTable$mdiff, 52.65),
#                            charge = c(adductsTable$charge, 1),
#                            n = c(adductsTable$n, 1),
#                            stringsAsFactors = F)

## ---- echo=T, eval=F-----------------------------------------------------
#  # The new adductsTable has to be also uploaded in the dbs list.
#  dbs <- assignDB()
#  dbs$adductsTable <- adductsTable
#  
#  idPCpos(MS1 = LipidMS::mix_pos_fullMS, MSMS1 = LipidMS::mix_pos_Ce20,
#  MSMS2 = LipidMS::mix_pos_Ce40, adducts = c("M+H", "M+Na", "M+X"),
#  dbs = dbs)

