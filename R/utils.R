#' Drop specified metabolites from data
#'
#' @param toDrop Metabolite to be removed
#' @param theData Data to remove metabolites from
#'
#' @returns `theData` without `toDrop` metabolites
#'
#' @importFrom collapse fsubset
#'
#' @export
#' @keywords internal
dropMetabolites <- function(toDrop, theData) {
  Metabolite <- NULL
  fsubset(theData, Metabolite %notin% toDrop)
}


#' Select the point geom based on interactivity
#'
#' @param interactive Whether the point should be interactive
#'
#' @returns `ggplot2::geom_point()` if interactive else `ggiraph::geom_point_interactive()`
#'
#' @keywords internal
getGeomPoint <- function(interactive) {
  if (interactive) {
    rlang::check_installed("ggiraph", reason = "ggiraph is required for interactive plots")
    return(ggiraph::geom_point_interactive)
  }
  ggplot2::geom_point
}
