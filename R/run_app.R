#' Run the interactive Shiny application for Bayesian transition priors
#'
#' @return No return value, called for side effects (launches a Shiny app).
#' @export
#'
#' @examples
#' \dontrun{
#' raretrans::run_app()
#' }
run_app <- function() {
  appDir <- system.file("shiny", "bayesian_transitions", package = "raretrans")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing raretrans.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
