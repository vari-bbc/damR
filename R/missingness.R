#' Get metabolites where proportion present is less than a set threshold
#'
#' @param qcData QC'd dataset
#' @param presence_threshold A threshold for non-missing data below which a
#' metabolite should be dropped
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


#' Remove missing data
#'
#' @param qcData input data
#' @param dropMetabolite Names of metabolites to drop
#' @returns Data with specified metabolites removed
#'
#' @export
RemoveMissingness <- function(qcData, dropMetabolite) {
  Metabolite <- NULL
  fsubset(qcData, Metabolite %notin% dropMetabolite)
}


#' Plot missingness heatmap
#'
#' Selects the necessary columns to make a wide form data frame. Makes the
#' presence/absence matrix with the Sample column and numeric columns and an
#' annotation data frame for annotation colors using the Subgroups.
#'
#' @param qcData Input data that has been QC'd
#' @param anno.pal RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation colors
#' @returns ComplexHeatmap::pheatmap
#'
#' @importFrom collapse fmutate fselect get_vars num_vars pivot
#'
#' @export
PlotMissingness <- function(qcData, anno.pal = c("Set1", "Set2")) {
  value <- NULL

  cols <- c("Sample", paste0("Subgroup", 1:2), "Metabolite", "value")
  data.bin.wide <- fselect(qcData, cols) |>
    fmutate(value = ifelse(is.na(value), 0, 1)) |>
    pivot(how = "wider", names = "Metabolite")
  metab.cols <- num_vars(data.bin.wide, return = "names")
  pmat <- fselect(data.bin.wide, c("Sample", metab.cols)) |>
    as.matrix(rownames = "Sample") |>
    t()

  anno <- get_vars(data.bin.wide, "Subgroup", regex = TRUE)
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
