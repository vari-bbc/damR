#' Compute PCA and variance explained
#'
#' @param fullData The full dataset including non-numeric data
#' @param numericColumns Column indices of numeric data that is put through `prcomp`
#' @returns A named list where `pca` is the PCA data column bound to `fullData`
#' and `varex` which contains the percent variance explained for each PC.
#'
#' @importFrom stats prcomp
#'
#' @export
GetPCA <- function(fullData, numericColumns) {
  pca <- prcomp(fullData[, numericColumns])
  pca.df <- cbind(as.data.frame(pca$x), fullData)
  pca.varex <- summary(pca)$importance[2, ]
  list(pca = pca.df, varex = pca.varex * 100)
}

#' Plot PCA
#'
#' PCA scatter plot with blanks and QC pool observations marked with a border
#'
#' @param pca The PCA and variance explained data from `GetPCA()`
#' @param x,y The principal component numbers to plot on the X and Y axes: numeric
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
PlotPCA <- function(pca, x = 1, y = 2, colorBy = "Group", pointSize = 3, pointAlpha = 0.8) {
  stopifnot(
    "PCA must contain PCA data and variance explained - use GetPCA()" = {
      names(pca) == c("pca", "varex")
    }
  )

  pd <- pca$pca
  xaxis <- paste0("PC", x)
  yaxis <- paste0("PC", y)
  pd$pool_blank <- grepl("pool|blank", pd$Group, ignore.case = TRUE)

  ggplot(
    pca$pca,
    aes(x = .data[[xaxis]], y = .data[[yaxis]], fill = .data[[colorBy]])
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
#' @param blankFoldThreshold Fold change threshold below which data are
#' considered blanks
#' @param metadataCols Column index of metadata
#'
#' @export
BlankFilter <- function(TheData, blankFoldThreshold = 5, metadataCols = 6) {
  metaboliteCols <- (metadataCols + 1):ncol(TheData)
  threshold <- log2(blankFoldThreshold)

  emptyTable <- data.frame(
    Metabolite = character(),
    Samples = character(),
    stringsAsFactors = FALSE
  )

  # Process blanks follow the same "PB" prefix convention QCData uses to drop controls
  isBlank <- substr(TheData$Group, 1, 2) == "PB" |
    substr(TheData$SampleID, 1, 2) == "PB"
  isBlank[is.na(isBlank)] <- FALSE

  if (!any(isBlank)) {
    print("Blank filtering requested but no PB samples were found; skipping.")
    return(list(TheData = TheData, blankTable = emptyTable, blankCount = 0))
  }

  # Working copy for the comparison only. NAs count as 0 here, matching the
  # reference implementation; they stay NA in TheData itself.
  e <- as.matrix(TheData[, metaboliteCols, drop = FALSE])
  storage.mode(e) <- "numeric"
  e[is.na(e)] <- 0

  # Compare each sample against the blanks run in its own batch
  sampleBatches <- as.character(TheData$Batch)
  sampleBatches[sampleBatches %in% c("NA", "")] <- NA
  batches <- sort(unique(sampleBatches[!is.na(sampleBatches)]))

  if (length(batches) > 0 && any(isBlank & !is.na(sampleBatches))) {
    for (batch in batches) {
      batchRows <- !is.na(sampleBatches) & sampleBatches == batch
      batchBlankRows <- batchRows & isBlank

      # Reference version errors here; in the app fall back to the pooled blanks
      if (sum(batchBlankRows) < 1) {
        print(paste0(
          "No PB sample found for batch ",
          batch,
          "; using all blanks instead."
        ))
        batchBlankRows <- isBlank
      }

      blanks <- apply(e[batchBlankRows, , drop = FALSE], 2, max)
      e[batchRows, ] <- sweep(e[batchRows, , drop = FALSE], 2, blanks, "-")
    }
  } else {
    blanks <- apply(e[isBlank, , drop = FALSE], 2, max)
    e <- sweep(e, 2, blanks, "-")
  }

  # Flag every sample value that did not clear its blank by the required fold change
  lowSignal <- e < threshold
  lowSignal[is.na(lowSignal)] <- FALSE
  lowSignal[isBlank, ] <- FALSE

  ## Apply to the real data ----
  beforeBlank <- as.matrix(TheData[, metaboliteCols, drop = FALSE])
  afterBlank <- beforeBlank
  afterBlank[lowSignal] <- NA_real_
  TheData[, metaboliteCols] <- afterBlank

  ## Table of what was removed ----
  # Report only values that actually went from a real number to NA, the way
  # blank_peptides() does in "Updated Blank Filter.R". Cells that were already NA
  # also fall below the threshold, since they count as 0 in the comparison above,
  # so indexing off lowSignal would pad the table with values that were never
  # there to remove.
  actuallyRemoved <- !is.na(beforeBlank) & is.na(afterBlank)
  removedIdx <- which(actuallyRemoved, arr.ind = TRUE)

  if (nrow(removedIdx) > 0) {
    removed <- data.frame(
      Metabolite = colnames(beforeBlank)[removedIdx[, "col"]],
      Sample = as.character(TheData$Sample)[removedIdx[, "row"]],
      stringsAsFactors = FALSE
    )

    # Base R rather than dplyr::summarise, since plyr is loaded after dplyr and masks it
    collapsed <- tapply(removed$Sample, removed$Metabolite, function(x) {
      paste(x, collapse = ", ")
    })
    blankTable <- data.frame(
      Metabolite = names(collapsed),
      Samples = as.character(collapsed),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    blankTable <- emptyTable
  }

  print(paste0(
    "Blank filtering set ",
    nrow(removedIdx),
    " values across ",
    nrow(blankTable),
    " metabolites to NA."
  ))

  return(list(
    TheData = TheData,
    blankTable = blankTable,
    blankCount = nrow(removedIdx)
  ))
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
  blankTable <- data.frame(
    Metabolite = character(),
    Samples = character(),
    stringsAsFactors = FALSE
  )
  blankCount <- 0

  if (isTRUE(blankFilter)) {
    blankList <- BlankFilter(TheData, blankFoldThreshold = blankFoldThreshold)
    TheData <- blankList$TheData
    blankTable <- blankList$blankTable
    blankCount <- blankList$blankCount
  }

  # Remove any groups and SampleIDs with TC, PB, or QC as first two letters
  TheData <- TheData[!substr(TheData$Group, 1, 2) %in% c("TC", "PB", "QC"), ]
  TheData <- TheData[!substr(TheData$SampleID, 1, 2) %in% c("TC", "PB", "QC"), ]

  # Trim metabolite numbers whitespace and force numeric
  TheData[, metaboliteCols] <- sapply(TheData[, metaboliteCols], function(x) {
    as.numeric(trimws(x))
  })

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
    blankCount = blankCount
  ))
}
