# sepByCE
#' Separate .mzXML files by CE
#'
#' Separation of .mzXML files from all-ions data by collision energy to work
#' with them separately.
#'
#' @param file path of the input .mzXML file
#' @param output a unique character value indicating the name of the output
#' files. The energy employed and .mzXML will be added automatically to each
#' file.
#'
#' @return As many .mzXML files as different collision energies employed.
#'
#' @details This function has been designed based on mzXML files obtained from
#' .d files (Agilent) using msConvert tool, in which we can find the collision
#' energy information. In addition to separate files by collision energies, this
#' function also changes the MS level of the high energy scans from 2 to 1
#' allowing their treatment (peak-picking for each collision energy, alignment,
#' i.e) with common software (xcms, mzMine2, enviPick, etc).
#'
#' @note Be careful with input and output arguments. For example, "file.mzXML"
#' would be the input argument and "file_sep" could be the output.
#'
#' @examples
#'
#' \dontrun{
#' sepByCE("input_file.mzXML", "output_file")
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maialba@alumni.uv.es>
sepByCE <- function(file, output){
  if (file.exists(file)){
    lines <- readLines(file)
    MS1 <- grep("          msLevel=\"1\"", lines)
    if (length(MS1) > 0){
      fullMS <- TRUE
    } else {
      fullMS <- FALSE
    }
    MS2 <- grep("          msLevel=\"2\"", lines)
    lCe <- grep("          collisionEnergy=.+\"", lines)
    if (length(MS2) > 0 || length(lCe) > 0){
      if (length(lCe) > 0){
        Ce <- grep("          collisionEnergy=.+\"", lines, value = T)
        Ce <- unlist(strsplit(Ce, "          collisionEnergy="))
        Ce <- unlist(strsplit(Ce[Ce != ""], "\""))
        Ce <- as.numeric(Ce[Ce != ""])
        CE <- unique(Ce)
        lines <- gsub("          msLevel=\"2\"", "          msLevel=\"1\"", lines)
        scans <- grep("<scan", lines)
        runs <- grep("</msRun>", lines)
        v <- vector()
        if(fullMS == TRUE){
          v <- append(v, c(CE[1:length(which(MS2 < MS1[1]))], "fullMS"))
          v <- append(v, CE[length(which(MS2 < MS1[1]))+1:length(CE)])
          v <- v[which(v != "NA")]
        } else {
          v <- as.character(CE)
        }
        for (i in 1:length(v)){
          pos <- seq(i,length(scans), length(v))
          lines2write <- c(1:(scans[1]-1))
          for (x in pos[1:length(pos)-1]){
            if (x == scans[length(scans)]){
              lines2write <- append(lines2write, c(scans[x]:c(runs-1)))
            } else {
              lines2write <- append(lines2write, c(scans[x]:(scans[x+1]-1)))
            }
          }
          lines2write <- append(lines2write, c(runs:length(lines)))
          write(lines[lines2write], file=paste(c(output, v[i], ".mzXML"),
            collapse=""))
        }
      } else {
        stop("No collision energies are detailed in the file!")
      }
    } else {
      stop("No scans from MS level 2 are present in the file!")
    }
  } else {
    stop("The input file does not exist!")
  }
}
