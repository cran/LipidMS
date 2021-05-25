# LipidMSapp()
#' LipidMS shiny app
#'
#' Interactive UI for LipidMS
#'
#' @examples
#' \dontrun{
#' library(LipidMS)
#' LipidMSapp()
#' }
#'
#' @author M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
LipidMSapp <- function() {
  if (interactive()){
    shiny::runApp(system.file('LipidMSapp', package='LipidMS'))
  }else {
    stop("Only available for interactive sessions")
  }
}
