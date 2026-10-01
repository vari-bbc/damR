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

  emptyTable <- data.frame(Metabolite = character(), Samples = character(), stringsAsFactors = FALSE)

  # Process blanks are matched two ways so either naming convention is picked up:
  # a leading "PB" (the prefix QCData uses to drop controls) or the word "blank"
  # anywhere in the group or sample ID, which is how the exports actually name them
  # ("TC_Blank-process_1" parses to SampleID "TC", Group "Blank_process", so the
  # PB prefix alone never matched it). Case insensitive, and "PB" is anchored so it
  # can't hit a group that merely contains those letters.
  blankPattern <- "^PB|blank"
  isBlank <- grepl(blankPattern, TheData$Group, ignore.case = TRUE) |
    grepl(blankPattern, TheData$SampleID, ignore.case = TRUE)
  isBlank[is.na(isBlank)] <- FALSE

  if (!any(isBlank)) {
    print("Blank filtering requested but no process blank samples were found; skipping.")
    return(list(
      TheData = TheData,
      blankTable = emptyTable,
      blankCount = 0,
      blanksFound = FALSE,
      perBatch = FALSE,
      fallbackBatches = character()
    ))
  }

  # QCData drops every TC/PB/QC row immediately after this function returns, so
  # those rows must be left alone: filtering them would otherwise put removals for
  # samples that never reach the analysis (the QC pools especially) into the table.
  controlPrefixes <- c("TC", "PB", "QC")
  isControl <- isBlank |
    substr(TheData$Group, 1, 2) %in% controlPrefixes |
    substr(TheData$SampleID, 1, 2) %in% controlPrefixes

  # Working copy for the comparison only. NAs count as 0 here, matching the
  # reference implementation; they stay NA in TheData itself.
  e <- as.matrix(TheData[, metaboliteCols, drop = FALSE])
  storage.mode(e) <- "numeric"
  e[is.na(e)] <- 0

  # Compare each sample against the blanks run in its own batch.
  # Leading zeros are stripped first: samples carry batch "01"/"02"/"03" while the
  # controls carry "1"/"2"/"3", so matching them as plain strings meant no batch
  # ever contained a blank and every batch quietly fell back to the pooled blanks.
  sampleBatches <- trimws(as.character(TheData$Batch))
  sampleBatches[sampleBatches %in% c("NA", "")] <- NA
  # Only normalise purely numeric batches, so text batch names are left as typed
  numericBatch <- !is.na(sampleBatches) & grepl("^[0-9]+$", sampleBatches)
  sampleBatches[numericBatch] <- as.character(as.numeric(sampleBatches[numericBatch]))
  batches <- sort(unique(sampleBatches[!is.na(sampleBatches)]))

  # Batches that had no blank of their own, so the pooled blanks stood in. Only
  # worth reporting for batches holding real samples; a batch of controls only
  # gets dropped regardless.
  fallbackBatches <- character()

  # FALSE when there is no usable batch information at all and every sample is
  # compared against the pooled blanks
  perBatch <- length(batches) > 0 && any(isBlank & !is.na(sampleBatches))

  if (perBatch) {
    for (batch in batches) {
      batchRows <- !is.na(sampleBatches) & sampleBatches == batch
      batchBlankRows <- batchRows & isBlank

      # Reference version errors here; in the app fall back to the pooled blanks
      if (sum(batchBlankRows) < 1) {
        print(paste0("No process blank found for batch ", batch, "; using all blanks instead."))
        if (any(batchRows & !isControl)) {
          fallbackBatches <- c(fallbackBatches, batch)
        }
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
  lowSignal[isControl, ] <- FALSE

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
    collapsed <- tapply(removed$Sample, removed$Metabolite, function(x) paste(x, collapse = ", "))
    blankTable <- data.frame(
      Metabolite = names(collapsed),
      Samples = as.character(collapsed),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    blankTable <- emptyTable
  }

  print(paste0("Blank filtering set ", nrow(removedIdx), " values across ", nrow(blankTable), " metabolites to NA."))

  return(list(
    TheData = TheData,
    blankTable = blankTable,
    blankCount = nrow(removedIdx),
    blanksFound = TRUE,
    perBatch = perBatch,
    fallbackBatches = fallbackBatches
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
