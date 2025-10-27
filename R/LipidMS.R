.onAttach <- function(lib, pkg){
  packageStartupMessage("
  LipidMS version 3.1.0: Lipid Annotation for LC-MS/MS DIA and DDA data analysis.
   --> Additional examples can be found at `github.com/maialba3/LipidMS/tree/main/inst/extdata`
   --> Batch processing is now available: batchdataProcessing() and annotatemsbatch().
   --> Shiny app for automatic processing and annotation through LipidMSapp().
   --> RT modelling to detect incorrect annotations and to propose new ones based on previous high confidence annotations based on fragmentation rules
     -> Interactive shiny app to check the proposed annotations based on RT criteria through manualfilterapp().
     -> Interactive shiny app to check new predicted annotations based on RT through manualpredictionapp(). ")
}
