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

#' Create the data and annotations for heatmaps
#'
#' @details
#' Used for plotting missingness and the significant metabolite heatmap.
#'
#' @param dt The long-form metabolite data from qc for missingness or imputed
#' data that's subset to significant metabolites.
#' @param anno.pal RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation colors
#'
#' @returns list with Data matrix, annotation data frame, and colors for the subgroups
#'
#' @keywords internal
getHeatmapInputs <- function(dt, anno.pal) {
  cols <- c("Sample", paste0("Subgroup", 1:2), "Metabolite", "value")
  data.bin.wide <- fselect(dt, cols) |>
    pivot(how = "wider", names = "Metabolite")
  metab.cols <- num_vars(data.bin.wide, return = "names")
  pmat <- fselect(data.bin.wide, c("Sample", metab.cols)) |>
    as.matrix(rownames = "Sample") |>
    t()
  anno <- get_vars(data.bin.wide, "Subgroup", regex = TRUE)
  anno.colors <- .getAnnoColors(anno, anno.pal)
  list(pmat, anno, anno.colors)
}

#' Helper to set up annotation colors for Sample groups in heatmap
#'
#' @param anno.df Annotation `data.frame` containing Subgroups as columns
#' @param anno.pal RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation colors
#' @returns Named list of named vectors of colors. List names are `anno.df`
#' column names and vector names are of the unique elements in `anno.df`.
#'
#' @examples
#' (anno.df <- data.frame(Subgroup1 = rep(letters[1:2], 3), Subgroup2 = letters[3:8]))
#' .getAnnoColors(anno.df, c("Set1", "Set2"))
#'
#' @export
#' @keywords internal
.getAnnoColors <- function(anno.df, anno.pal) {
  names(anno.pal) <- colnames(anno.df)
  sapply(
    names(anno.df),
    function(x) {
      counts <- length(unique(anno.df[[x]]))
      pal <- RColorBrewer::brewer.pal(counts, anno.pal[x])
      names(pal) <- unique(anno.df[[x]])
      pal[seq(counts)]
    },
    simplify = FALSE
  )
}
