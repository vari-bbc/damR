#' Plot heatmap of significantly differentially abundant metabolites
#'
#' @details
#' Additional plot aesthetics can be controlled through arguments passed to
#' `ComplexHeatmap::pheatmap()` though `...`. However, `mat`, `annotation_col`,
#' and `annotation_colors` are set by `damR` so modifying those arguments may
#' fail. Use `anno.pal` to modify the annotation colors.
#'
#' @param diffMetabolites The differentially abundant metabolite data from
#' `DiffMetabolites()` used to select metabolites for the heatmap
#' @param imputedData Imputed metabolite data from `impute.QRILC()` used to
#' produce the heatmap
#' @param fcThreshold The fold change to filter metabolites as 'significant'
#' @param alpha Significance level to filter the differential abundance p-value
#' @param scale Whether to scale the imputed data matrix to mean 0 and standard
#' deviation of 1
#' @param anno.pal RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation colors
#' @inheritDotParams ComplexHeatmap::pheatmap -mat -annotation_col -annotation_colors
#'
#' @returns ComplexHeatmap::pheatmap
#'
#' @importFrom collapse fsubset %=% .c
#'
#' @export
PlotSigMetabolites <- function(
  diffMetabolites,
  imputedData,
  fcThreshold = 2,
  alpha = 0.05,
  scale = FALSE,
  anno.pal = c("Set1", "Set2"),
  ...
) {
  logFC <- adj.P.Val <- Metabolite <- pmat <- anno <- anno.colors <- NULL

  sigMetabolites <- fsubset(diffMetabolites, logFC >= fcThreshold & adj.P.Val < alpha) |>
    _[["Compound"]] |>
    funique()
  sigImputed <- fsubset(imputedData, Metabolite %in% sigMetabolites)
  .c(pmat, anno, anno.colors) %=% getHeatmapInputs(sigImputed, anno.pal)
  if (scale) {
    pmat <- scale(pmat)
  }
  .plotSigHeatmap(
    pmat,
    scaled = scale,
    anno,
    anno.colors,
    ...
  )
}

#' Helper to plot missingness heatmap
#'
#' @param mat Matrix of imputed data subset to differentially abundant metabolites
#' @param scaled Whether the input matrix was scaled
#' @param anno.df Annotation data.frame with the sample groups
#' @param anno.pal List named with `anno.df` column names containing named
#' color vectors where the color vector names correspond to the elements of
#' anno. See `.getAnnoColors()`
#' @inheritDotParams ComplexHeatmap::pheatmap -mat -annotation_col -annotation_colors
#'
#' @returns ComplexHeatmap
#'
#' @export
#' @keywords internal
.plotSigHeatmap <- function(
  mat,
  scaled,
  anno.df,
  anno.pal,
  ...
) {
  htitle <- ifelse(scaled, expression("Scaled" ~ log[2] ~ "Abundance"), expression(log[2] ~ "Abundance"))
  ComplexHeatmap::pheatmap(
    mat,
    annotation_col = anno.df,
    heatmap_legend_param = list(title = htitle, title_position = "leftcenter-rot"),
    annotation_colors = anno.pal,
    ...
  )
}
