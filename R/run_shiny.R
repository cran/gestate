#'Load Shiny for Gestate
#' Loads the Shiny interactive GUI for gestate
#' @examples \donttest{run_gestate()}
#' @export
#' @import shiny
#' @import shinythemes
run_gestate <- function() {
  appDir <- system.file("Shiny", package = "gestate")
  if (appDir == "") {
    stop("Could not find Shiny installation directory. Check that the Shiny folder exists within the gestate installation folder. To fix, try re-installing `gestate`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
