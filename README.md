# LipidMS v3.1
Lipid Annotation for LC-MS/MS DIA and DDA data analysis. New features:

  - Batch processing: peak-picking, grouping and alignment wrapped in batchdataProcessing(). Lipid annotation for msbatch objects simplified with annotatemsbatch().
  - New lipid classes: plasmanyl and plasmenyl PC and PE, acylceramides and ceramides phosphate.
  - GUI through shiny app running LipidMSapp().
  - Improved graphical outputs for lipid annotation.
  - RT modelling to detect incorrect annotations and to propose new ones based on previous high confidence annotations identified by fragmentation rules.
    - Interactive shiny app to check the proposed annotations based on RT criteria through manualfilterapp().
    - Interactive shiny app to check new predicted annotations based on RT through manualpredictionapp(). 


Citation:

If you use this software in your research, please cite:

  - Alcoriza-Balaguer MI., et al. (2019) LipidMS: An R Package for Lipid Annotation in Untargeted Liquid Chromatography-Data Independent Acquisition-Mass Spectrometry Lipidomics. Anal Chem, 2019, 91(1), 836-845. doi:10.1021/acs.analchem.8b03409.
  - LipidMS 3.0: an R-package and a web-based tool for LC-MS/MS data processing and lipid annotation. Bioinformatics, 2022. doi.org/10.1093/bioinformatics/btac581


References:

  - Peak-picking algorithm has been imported from enviPick R-package (Martin Loos) (https://github.com/blosloos/enviPick)
