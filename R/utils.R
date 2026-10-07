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
