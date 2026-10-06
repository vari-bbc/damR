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


#' Blank filter
#'
#' @details
#' Adapted from the proteomics blank_filter_by_batch() in "Updated Blank
#' Filter.R". Same logic, transposed layout: here samples are rows and
#' metabolites are columns. For each metabolite, a sample value must exceed the
#' loudest process blank in its own batch by at least blankFoldThreshold fold;
#' values that don't are set to NA. The blanks are only used to decide what to
#' drop -- surviving values keep their original intensities, they are not
#' blank-subtracted. Expects TheData metabolite columns to already be log2
#' transformed (as they are by the time QCData runs), which makes the fold
#' change a simple subtraction.
#'
#' @param TheData input data
#' @param blankFilter Whether to run blank filtering
#' @param blankFoldThreshold Fold change threshold below which data are
#' considered blanks
#'
#' @export
QCData <- function(TheData, blankFilter = FALSE, blankFoldThreshold = 5) {
  metaboliteCols <- 7:ncol(TheData)

  # ensure all col names are acceptable col names
  colnames(TheData) <- make.names(colnames(TheData))

  ## Blank filtration ----
  # Must run before the PB rows are dropped below, since it needs the blanks
  blankTable <- data.frame(Metabolite = character(), Samples = character(), stringsAsFactors = FALSE)
  blankCount <- 0
  # Tracks whether any blanks were actually found, so the app can tell "the filter
  # ran and nothing failed the threshold" apart from "the filter never ran"
  blanksFound <- FALSE
  perBatch <- FALSE
  fallbackBatches <- character()

  if (isTRUE(blankFilter)) {
    blankList <- BlankFilter(TheData, blankFoldThreshold = blankFoldThreshold)
    TheData <- blankList$TheData
    blankTable <- blankList$blankTable
    blankCount <- blankList$blankCount
    blanksFound <- isTRUE(blankList$blanksFound)
    perBatch <- isTRUE(blankList$perBatch)
    fallbackBatches <- blankList$fallbackBatches
  }

  # Remove any groups and SampleIDs with TC, PB, or QC as first two letters
  TheData <- TheData[!substr(TheData$Group, 1, 2) %in% c("TC", "PB", "QC"), ]
  TheData <- TheData[!substr(TheData$SampleID, 1, 2) %in% c("TC", "PB", "QC"), ]

  # Trim metabolite numbers whitespace and force numeric
  TheData[, metaboliteCols] <- sapply(TheData[, metaboliteCols], function(x) as.numeric(trimws(x)))

  # Find Zero Variance function
  zeroVar <- function(x) {
    (I(sd(x, na.rm = T) < 0.01 | is.na(sd(x, na.rm = T))))
  }

  # Get colnames of metabolites with 0 variance
  rmv <- apply(TheData[, metaboliteCols], 2, function(x) zeroVar(x))
  rmv <- rmv[rmv]
  rmvCols <- names(rmv)

  print("0 Variance Metabolites:")
  print(rmvCols)

  # remove them
  TheData <- TheData[, !colnames(TheData) %in% rmvCols]

  return(list(
    TheData = TheData,
    blankTable = blankTable,
    blankCount = blankCount,
    blanksFound = blanksFound,
    perBatch = perBatch,
    fallbackBatches = fallbackBatches
  ))
}
