#' Plot volcano for a single contrast
#'
#' @param resultTable The `topTable` results from `pAdjustData()`
#' @param contrastGroups Two contrast group names
#' @param contrastTitle Optional title for the Volcano plot
#' @param alphaLevel Significance cutoff level
#' @param cap Absolute log fold-change cap to cap the `logFC` x axis so that `-cap < logFC < cap`.
#' @returns Volcano ggplot
#'
#' @export
PlotVolcano <- function(
  resultTable,
  contrastGroups,
  contrastTitle = NULL,
  alphaLevel = 0.05,
  cap = 10
) {
  pd <- VolcanoData(resultTable, contrastGroups, contrastTitle, alphaLevel, cap)
  VolcanoPlotter(pd, alphaLevel) +
    labs(
      title = contrastTitle,
      caption = sprintf("all logFC capped at %1$.2f/-%1$.2f for plotting purposes", cap)
    )
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
VolcanoData <- function(vpd, contrastGroups, contrastTitle, alphaLevel = 0.05, cap = 10) {
  contrastStr <- MakeTheContrasts(contrastGroups, contrastTitle)
  fsubset(vpd, Contrast == contrastStr) |>
    fmutate(
      logFC = pmax(pmin(logFC, cap), -cap), # cap between -cap and cap
      Direction = fcase(
        logFC > 0 & adj.P.Val < alphaLevel , "Up regulated"   ,
        logFC < 0 & adj.P.Val < alphaLevel , "Down regulated" ,
        default = "NS"
      )
    )
}

#' Volcano plot
#'
#' @param pd Plot data produced by `VolcanoData()`
#' @param alphaLevel Significance cutoff level
#' @returns Volcano ggplot
#'
#' @importFrom ggiraph geom_point_interactive girafe girafe_options opts_selection
#' @importFrom ggplot2 aes facet_wrap geom_hline ggplot ggtitle labs scale_color_manual theme theme_minimal ylab
#' @importFrom rlang .data
#'
#' @export
VolcanoPlotter <- function(pd, alphaLevel) {
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
    ggiraph::geom_point_interactive(size = 2) +
    geom_hline(
      yintercept = -log10(alphaLevel),
      linetype = "dotted",
      color = "red",
      linewidth = 0.5
    ) +
    scale_color_manual(
      values = c(
        "Up regulated" = "#276419",
        "Down regulated" = "#8E0152",
        "NS" = "grey"
      )
    ) +
    labs(
      y = expression(-log[10] ~ "(FDR p-value)"),
      subtitle = pd$Contrast[1]
    ) +
    theme_minimal(7) +
    theme(legend.position = "top")
}

#' Plot volcanoes for all input contrasts
#'
#' @inheritParams PlotVolcano
#'
#' @returns A list of plots
#'
#' @importFrom ggiraph girafe_css opts_selection opts_sizing
#'
#' @export
PlotAllVolcanoes <- function(
  resultTable,
  contrastGroups,
  contrastTitle,
  alphaLevel = 0.05,
  cap = 10
) {
  compareGroups <- combn(contrastGroups, 2, simplify = FALSE)
  lapply(compareGroups, \(x) PlotVolcano(resultTable, x, contrastTitle, alphaLevel, cap))
}
