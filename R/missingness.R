#' Get metabolites with missingness proportion below a set threshold
#'
#' @param qcData QC'd dataset
#' @param presence_threshold A threshold for non-missing data below which a
#' metabolite should be dropped
#' @returns Vector of metabolite names below the threshold
#'
#' @export
IdentifyMissing <- function(qcData, presence_threshold = 0.75) {
  metaboliteCols <- 7:ncol(qcData)

  metab.mat <- as.matrix(qcData[, metaboliteCols])
  metab.mat[metab.mat == 0] <- NA # Turn 0s back into NA -?
  pa.mat <- ifelse(is.na(metab.mat), 0, 1)

  pPresent <- apply(pa.mat, 2, function(x) sum(x) / length(x))
  names(pPresent[pPresent < presence_threshold])
}


#' Remove missing data
#'
#' @param qcData input data
#' @param dropCols Column names of metabolites to drop
#' @returns A list of data with missing values removed, and missingness heatmap and table
#'
#' @export
RemoveMissingness <- function(qcData, dropCols) {
  qcData[, !colnames(qcData) %in% dropCols]
}


#' Plot missingness heatmap
#'
#' @param qcData Input data that has been QC'd
#' @param anno.pal RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation colors
#'
#' @returns ComplexHeatmap::pheatmap
#' @export
PlotMissingness <- function(qcData, anno.pal = c("Set1", "Set2")) {
  pmat <- qcData[, 7:ncol(qcData)] |> as.matrix() |> t()
  pmat[!is.na(pmat)] <- 1
  pmat[is.na(pmat)] <- 0
  colnames(pmat) <- qcData$Sample

  anno <- qcData[, c("Subgroup1", "Subgroup2")]
  rownames(anno) <- qcData$Sample

  anno.colors <- .getAnnoColors(anno, anno.pal)
  .plotMissingHeatmap(pmat, anno, anno.colors)
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
    simplify = F
  )
}

#' Helper to plot missingness heatmap
#'
#' @param bin.mat Binary matrix where 0 is missing and 1 is present
#' @param anno.df Annotation data.frame with the sample groups
#' @param anno.colors List named with `anno.df` column names containing named
#' color vectors where the color vector names correspond to the elements of
#' anno. See `.getAnnoColors()`
.plotMissingHeatmap <- function(bin.mat, anno.df, anno.colors) {
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
    show_rownames = FALSE,
    annotation_colors = anno.colors
  )
}
