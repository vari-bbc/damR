#' Plot volcano for a single contrast
#'
#' @param resultTable The `topTable` results from `pAdjustData()`
#' @param contrastGroups Two contrast group names
#' @param contrastTitle Optional title for the Volcano plot
#' @param alphaLevel Significance cutoff level
#' @param cap Absolute log fold-change cap to cap the `logFC` x axis so that `-cap < logFC < cap`.
#' @param ... Geom aesthetics passed to `VolcanoPlotter()`.
#' @returns Volcano ggplot
#'
#' @export
PlotVolcano <- function(
  resultTable,
  contrastGroups,
  contrastTitle = NULL,
  alphaLevel = 0.05,
  cap = NULL,
  ...
) {
  stopifnot("`contrastGroups` accepts only two groups from the data" = length(contrastGroups) == 2)
  pd <- VolcanoData(resultTable, contrastGroups, contrastTitle, alphaLevel, cap)
  VolcanoPlotter(pd, alphaLevel, ...) +
    labs(
      title = contrastTitle,
      caption = ifelse(is.null(cap), "", sprintf("all logFC capped at %1$.2f/-%1$.2f for plotting purposes", cap))
    ) +
    theme(plot.caption = element_text(size = rel(0.5)))
}

#' Produce volcano plot data for a single contrast
#'
#' @inheritParams PlotVolcano
#'
#' @returns data.table ready for plotting with `VolcanoPlotter()`
#'
#' @importFrom collapse fsubset fmutate
#' @importFrom data.table fcase
#'
#' @export
VolcanoData <- function(resultTable, contrastGroups, contrastTitle, alphaLevel = 0.05, cap = NULL) {
  Contrast <- logFC <- adj.P.Val <- NULL
  contrastStr <- MakeTheContrasts(contrastGroups, contrastTitle)
  vd <- fsubset(resultTable, Contrast == contrastStr) |>
    fmutate(
      Direction = fcase(
        logFC > 0 & adj.P.Val < alphaLevel , "Up regulated"   ,
        logFC < 0 & adj.P.Val < alphaLevel , "Down regulated" ,
        default = "NS"
      )
    )
  if (is.null(cap)) {
    return(vd)
  }
  vd <- fmutate(vd, logFC = pmax(pmin(logFC, cap), -cap))
}

#' Volcano plot
#'
#' @param pd Plot data produced by `VolcanoData()`
#' @param alphaLevel Significance cutoff level
#' @param pointColors Named vector of colors for `"Up regulated"`, `"Down
#' regulated"`, and `"NS"` points
#' @param pointSize Size of the `ggplot2::geom_point()`
#' @param pointAlpha The transparency of the [`ggplot2::geom_point()`]
#' @param alphaLineColor,alphaLineType,alphaLineWidth Color, type and width of
#' the line used to mark the significance cutoff. See `?ggplot2::geom_hline()`.
#' @param interactive Whether to plot as interactive with `ggiraph::geom_point_interactive()`
#' @returns Volcano ggplot
#'
#' @importFrom ggplot2 ggplot aes facet_wrap geom_hline geom_point ggtitle labs scale_color_manual theme theme_minimal ylab
#' @importFrom rlang .data
#'
#' @export
#' @keywords internal
VolcanoPlotter <- function(
  pd,
  alphaLevel,
  pointColors = c(
    "Up regulated" = "#276419",
    "Down regulated" = "#8E0152",
    "NS" = "grey"
  ),
  pointSize = 2,
  pointAlpha = 0.7,
  alphaLineColor = "red",
  alphaLineType = "dotted",
  alphaLineWidth = 0.5,
  interactive = FALSE
) {
  point_aes <- getGeomPoint(interactive)

  ggplot(
    pd,
    aes(
      x = .data$logFC,
      y = -log10(.data$adj.P.Val),
      color = .data$Direction,
      data_id = .data$Compound,
      tooltip = .data$Compound,
      label = .data$Compound
    )
  ) +
    point_aes(size = pointSize, alpha = pointAlpha) +
    geom_hline(
      yintercept = -log10(alphaLevel),
      linetype = alphaLineType,
      color = alphaLineColor,
      linewidth = alphaLineWidth
    ) +
    scale_color_manual(values = pointColors) +
    labs(y = expression(-log[10] ~ "(FDR p-value)"), subtitle = pd$Contrast[1]) +
    theme(legend.position = "top")
}

#' Plot volcanoes for all input contrasts
#'
#' @inheritParams PlotVolcano
#'
#' @returns A list of plots
#'
#' @importFrom stats setNames
#'
#' @export
PlotAllVolcanoes <- function(
  resultTable,
  contrastGroups,
  contrastTitle,
  alphaLevel = 0.05,
  cap = NULL,
  ...
) {
  compareGroups <- combn(contrastGroups, 2, simplify = FALSE)
  plot_names <- vapply(compareGroups, \(x) paste(x, collapse = " - "), FUN.VALUE = "character")
  plots <- lapply(compareGroups, \(x) PlotVolcano(resultTable, x, contrastTitle, alphaLevel, cap, ...))
  setNames(plots, plot_names)
}
