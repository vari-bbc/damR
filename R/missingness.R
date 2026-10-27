#' Remove missing data
#'
#' @param qcData QC'd dataset
#' @param presence_threshold A threshold for non-missing data below which a
#' metabolite should be dropped
#'
#' @returns Data with specified metabolites removed
#'
#' @export
RemoveMissing <- function(qcData, presence_threshold = 0.75) {
  toDrop <- IdentifyMissing(qcData, presence_threshold)
  dropMetabolites(toDrop, qcData)
}

#' Get metabolites where proportion present is less than a set threshold
#'
#' @inheritParams RemoveMissing
#'
#' @returns Vector of metabolite names below the threshold
#'
#' @importFrom collapse fgroup_by fmean fmutate fsubset fsummarise fungroup recode_num
#'
#' @export
IdentifyMissing <- function(qcData, presence_threshold = 0.75) {
  value <- pPresent <- NULL
  noZero <- recode_num(qcData, `0` = NA)
  prop.present <- fgroup_by(noZero, "Metabolite") |>
    fsummarise("pPresent" = 1 - fmean(is.na(value)))
  fsubset(prop.present, pPresent < presence_threshold)[["Metabolite"]]
}

#' Plot missingness heatmap
#'
#' Selects the necessary columns to make a wide form data frame. Makes the
#' presence/absence matrix with the Sample column and numeric columns and an
#' annotation data frame for annotation colors using the Subgroups.
#'
#' @param qcData Input data that has been QC'd
#' @param anno.pal RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation colors
#' @param showRownames Whether to show metabolite names as rownames in the heatmap
#' @returns ComplexHeatmap::pheatmap
#'
#' @importFrom collapse fmutate fselect get_vars num_vars pivot
#'
#' @export
PlotMissingness <- function(qcData, anno.pal = c("Set1", "Set2"), showRownames = TRUE) {
  value <- pmat <- anno <- anno.colors <- NULL

  data.bin <- fmutate(qcData, value = ifelse(is.na(value), 0, 1))
  .c(pmat, anno, anno.colors) %=% getHeatmapInputs(data.bin, anno.pal)
  .plotMissingHeatmap(pmat, anno, anno.colors, showRownames)
}

#' Helper to plot missingness heatmap
#'
#' @param bin.mat Binary matrix where 0 is missing and 1 is present
#' @param anno.df Annotation data.frame with the sample groups
#' @param anno.colors List named with `anno.df` column names containing named
#' color vectors where the color vector names correspond to the elements of
#' anno. See `.getAnnoColors()`
#' @param showRownames Whether to show metabolite names as rownames in the heatmap
#'
#' @returns ComplexHeatmap
#'
#' @export
#' @keywords internal
.plotMissingHeatmap <- function(bin.mat, anno.df, anno.colors, showRownames) {
  ComplexHeatmap::pheatmap(
    bin.mat,
    color = c("grey70", "green3"),
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    treeheight_row = 0,
    annotation_col = anno.df,
    fontsize_row = 2,
    fontsize_col = 8,
    legend_breaks = c(0, 1),
    legend_labels = c("Missing", "Present"),
    heatmap_legend_param = list(
      title = "",
      at = c(0, 1),
      color_bar = "discrete"
    ),
    show_rownames = showRownames,
    annotation_colors = anno.colors
  )
}
