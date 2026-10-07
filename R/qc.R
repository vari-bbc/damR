#' Run QC on input data
#'
#' Removes:
#'   - data points below a fold change threshold compared to blank samples. See `BlankFilter()`.
#'   - metabolites with zero variance
#'   - Blank samples (after blank filtering)
#'   - Control samples
#'
#' @param inputData Input metabolite data
#' @param blankFilter Whether to run blank filtering
#' @param blankFoldThreshold The blank fold change threshold. The difference
#' between every metabolite and its batch's blank must be at least this value
#' or higher to remain in the dataset.
#' @param blankPattern String pattern to denote blank samples
#' @param controlPattern String pattern to denote control samples like QC
#'
#' @importFrom collapse fgroup_by fsubset fsummarise
#' @importFrom stats sd
#'
#' @export
QCData <- function(
  inputData,
  blankFilter = TRUE,
  blankFoldThreshold = 5,
  blankPattern = "^PB|blank",
  controlPattern = "^PB|^TC|^QC"
) {
  if (blankFilter) {
    inputData <- BlankFilter(inputData, blankFoldThreshold, blankPattern)
  }

  zeroVarMetabolites <- getZeroVar(inputData)
  zeroVarDropped <- dropMetabolites(zeroVarMetabolites, inputData)

  Sample <- NULL
  fsubset(inputData, !(Sample %like% controlPattern))
}

#' Identify metabolites with zero variance, including any that are all NA
#'
#' @param inputData Input metabolite data
#'
#' @importFrom collapse fgroup_by fsubset fsummarise
#' @importFrom stats sd
#'
#' @export
#' @keywords internal
getZeroVar <- function(inputData) {
  value <- zv <- NULL
  fgroup_by(inputData, "Metabolite") |>
    fsummarise(zv = sd(value, na.rm = TRUE) < 0.01) |>
    fsubset(zv == TRUE | is.na(zv)) |>
    _[["Metabolite"]]
}

#' Produce PCA plots for QC
#'
#' @param TheData input data
#' @param xAxis input data
#' @param yAxis input data
#' @param dimQCBlank input data
#' @param dimAlpha input data
#' @param allGroups input data
#' @returns plotly PCA plot
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics layout
#' @importFrom plotly add_trace plot_ly
#' @importFrom stats lm prcomp qnorm quantile sd setNames
#' @importFrom utils read.csv
#'
#' @export
#' @keywords internal
QCPCAPlots <- function(
  TheData,
  xAxis,
  yAxis,
  dimQCBlank = FALSE,
  dimAlpha = 0.3,
  allGroups = NULL
) {
  metaboliteCols <- 7:ncol(TheData)

  TheData[is.na(TheData)] = 0

  pc <- prcomp(TheData[, metaboliteCols])
  pc.df <- as.data.frame(pc$x)
  pc.df$Group <- TheData$Group

  pc.df$XAxis <- pc.df[, xAxis]
  pc.df$YAxis <- pc.df[, yAxis]

  pc.df <- pc.df
  #
  # PreQCPCAPlot <- ggplot(pc.df, aes(x = XAxis, y = YAxis, color = Group )) +
  #   geom_point(size = 0.5) +
  #   theme_minimal() + scale_color_viridis_d() +
  #   theme(legend.text = element_text(size = 2),
  #         legend.title = element_text(size = 4),
  #         legend.key.size = unit(0.25, "pt"),
  #         axis.text = element_text(size = 2),
  #         axis.title = element_text(size = 4)) +
  #   labs(x = xAxis, y = yAxis)

  pc.df$Group <- as.character(pc.df$Group)

  if (is.null(allGroups)) {
    allGroups <- sort(unique(pc.df$Group))
  } else {
    allGroups <- unique(as.character(allGroups))
  }

  pc.df$Group <- factor(pc.df$Group, levels = allGroups)

  groupLevels <- levels(pc.df$Group)
  cols <- viridis::viridis(length(groupLevels))
  groupColorMap <- setNames(cols, groupLevels)

  PreQCPCAPlot <- plot_ly(type = "scatter", mode = "markers")

  if (dimQCBlank) {
    qcBlankMask <- grepl(
      "pool|blank",
      as.character(pc.df$Group),
      ignore.case = TRUE
    )
    qcBlankDf <- pc.df[qcBlankMask, , drop = FALSE]

    if (nrow(qcBlankDf) > 0) {
      PreQCPCAPlot <- PreQCPCAPlot |>
        add_trace(
          x = qcBlankDf[[xAxis]],
          y = qcBlankDf[[yAxis]],
          marker = list(
            size = 11,
            color = "black"
          ),
          showlegend = FALSE,
          hoverinfo = "skip"
        )
    }
  }

  for (currentGroup in groupLevels) {
    groupDf <- pc.df[pc.df$Group == currentGroup, , drop = FALSE]

    if (nrow(groupDf) == 0) {
      next
    }

    PreQCPCAPlot <- PreQCPCAPlot |>
      add_trace(
        x = groupDf[[xAxis]],
        y = groupDf[[yAxis]],
        name = currentGroup,
        marker = list(
          size = 8,
          color = unname(groupColorMap[[currentGroup]])
        ),
        showlegend = TRUE
      )
  }

  PreQCPCAPlot <- PreQCPCAPlot |>
    plotly::layout(
      xaxis = list(
        title = xAxis,
        tickfont = list(size = 10),
        titlefont = list(size = 12)
      ),
      yaxis = list(
        title = yAxis,
        tickfont = list(size = 10),
        titlefont = list(size = 12)
      ),
      legend = list(
        font = list(size = 10)
      ),
      template = "plotly_white"
    )

  # PreQCPCAPlot <- p1 + p2 + p3 + plot_layout(ncol=3)

  return(PreQCPCAPlot)
}
