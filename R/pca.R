#' Compute PCA and variance explained
#'
#' @param fullData The full dataset including non-numeric data
#' @returns A named list where `pca` is the PCA data column bound to `fullData`
#' and `varex` which contains the percent variance explained for each PC.
#'
#' @importFrom data.table as.data.table
#' @importFrom collapse fselect num_vars pivot replace_na
#' @importFrom stats prcomp
#'
#' @export
GetPCA <- function(fullData) {
  data.wide <- fullData[, c("Sample", "Group", "Metabolite", "value")] |>
    replace_na(value = 0) |>
    pivot(how = "wider", ids = c("Sample", "Group"), names = "Metabolite")

  metab.cols <- num_vars(data.wide, return = "names")
  mat <- fselect(data.wide, c("Sample", metab.cols)) |>
    as.matrix(rownames = "Sample")

  pca <- prcomp(mat)
  pca.df <- cbind(
    data.wide[, c("Group")],
    as.data.table(pca$x, keep.rownames = "Sample")
  )
  pca.varex <- summary(pca)$importance[2, ]
  list(pca = pca.df, varex = pca.varex * 100)
}

#' Plot PCA
#'
#' PCA scatter plot with blanks and QC pool observations marked with a border
#'
#' @param pca The PCA and variance explained data from `GetPCA()`
#' @param x,y The principal components (prefixed with "PC") to plot on the X and Y axes
#' @param colorBy The column by which to color the points (using `aes(fill = ...)`)
#' @param pointSize The size of the [`ggplot2::geom_point()`]
#' @param pointAlpha The transparency of the [`ggplot2::geom_point()`]
#' @returns A named list where `pca` is the PCA data column bound to `fullData`
#' and `varex` which contains the percent variance explained for each PC.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs
#' @importFrom rlang .data
#'
#' @export
PlotPCA <- function(
  pca,
  x = "PC1",
  y = "PC2",
  colorBy = "Group",
  pointSize = 3,
  pointAlpha = 0.8
) {
  stopifnot(
    "PCA must contain PCA data and variance explained - use GetPCA()" = {
      names(pca) == c("pca", "varex")
    }
  )

  pd <- pca$pca
  xaxis <- x
  yaxis <- y
  pd$pool_blank <- grepl("pool|blank", pd$Group, ignore.case = TRUE)

  ggplot(
    pca$pca,
    aes(
      x = .data[[xaxis]],
      y = .data[[yaxis]],
      fill = .data[[colorBy]]
    )
  ) +
    geom_point(
      shape = 21,
      size = pointSize,
      alpha = pointAlpha,
      stroke = ifelse(pd$pool_blank == TRUE, 1, 0)
    ) +
    labs(
      x = sprintf("PCA (%.2f%%)", pca$varex[x]),
      y = sprintf("PCA (%.2f%%)", pca$varex[y])
    )
}
