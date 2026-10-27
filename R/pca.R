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
#' @param interactive Whether to plot as interactive with `ggiraph::geom_point_interactive()`
#'
#' @returns A named list where `pca` is the PCA data column bound to `fullData`
#' and `varex` which contains the percent variance explained for each PC.
#'
#' @importFrom collapse fmutate
#' @importFrom ggplot2 ggplot aes geom_point labs guides
#' @importFrom rlang .data
#'
#' @export
PlotPCA <- function(
  pca,
  x = "PC1",
  y = "PC2",
  colorBy = "Group",
  pointSize = 3,
  pointAlpha = 0.8,
  interactive = FALSE
) {
  stopifnot(
    "PCA must contain PCA data and variance explained - use GetPCA()" = {
      names(pca) == c("pca", "varex")
    }
  )

  point_aes <- getGeomPoint(interactive)
  Group <- pool_blank <- NULL
  pd <- pca$pca
  xaxis <- x
  yaxis <- y
  is_pool_blank <- grepl("pool|blank", pd$Group, ignore.case = TRUE)
  pd <- fmutate(pd, "pool_blank" = ifelse(is_pool_blank, Group, "Experimental sample"))
  if (any(is_pool_blank)) {
    pool_blank_names <- pca$pca[is_pool_blank, funique(Group)]
    pd <- fmutate(pd, pool_blank = factor(pool_blank, levels = c("Experimental sample", pool_blank_names)))
  }

  ggplot(
    pd[!is_pool_blank],
    aes(
      x = .data[[xaxis]],
      y = .data[[yaxis]],
      color = .data[[colorBy]],
      shape = .data$pool_blank,
      data_id = .data$Group,
      tooltip = .data$Group,
      label = .data$Group
    )
  ) +
    point_aes(
      size = pointSize,
      alpha = pointAlpha
    ) +
    point_aes(
      data = pd[is_pool_blank],
      aes(
        color = NULL
      ),
      size = pointSize,
      alpha = pointAlpha
    ) +
    labs(
      x = sprintf("%s (%.2f%%)", x, pca$varex[x]),
      y = sprintf("%s (%.2f%%)", y, pca$varex[y]),
      shape = "Sample type"
    ) +
    guides(shape = ifelse(!any(is_pool_blank), "none", "legend"))
}
