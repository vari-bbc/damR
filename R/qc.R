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


#' Remove missing data
#'
#' @param TheData input data
#' @returns A list of data with missing values removed, and missingness heatmap and table
#'
#' @export
RemoveMissingness <- function(TheData) {
  TheData <- TheData
  metaboliteCols <- 7:ncol(TheData)

  # Turn 0s back into NA
  TheData[TheData == 0] <- NA

  # Make seperate dataframe for heatmap
  miss.df <- TheData

  # Set rownames to sample
  rownames(miss.df) <- miss.df$sample

  # Subset to metabolite cols
  miss.df <- miss.df[, metaboliteCols]

  # If missing set to one, if not set to 0
  miss.df[!is.na(miss.df)] <- 1
  miss.df[is.na(miss.df)] <- 0
  rownames(miss.df) <- TheData$Sample

  # Get percent missing
  missing <- apply(miss.df, 2, function(x) sum(x) / length(x))

  # Get subgroups as anno
  anno <- TheData[, colnames(TheData) %in% c("Subgroup1", "Subgroup2")]
  rownames(anno) <- TheData$Sample

  # Get colors for annots
  mycolors <- list()
  pals <- list(
    RColorBrewer::brewer.pal(n = 3, "Set1"),
    RColorBrewer::brewer.pal(n = 5, "Set2")
  )

  for (i in 1:ncol(anno)) {
    newCols <- colorRampPalette(pals[[i]])
    mycols <- newCols(length(unique(anno[, i])))
    names(mycols) <- unique(anno[, i])
    mycolors[[paste0("Subgroup", i)]] = mycols
  }

  # Display matrix uses Present = 0 and Missing = 1 for the legend/key.
  miss.plot.df <- ifelse(miss.df == 1, 0, 1)

  ## Missingness Heatmap ----
  # Make heatmap
  if (length(which(missing < 0.9)) > 0) {
    missingnessHeatmap <- ComplexHeatmap::pheatmap(
      t(miss.plot.df[, which(missing < 0.9)]),
      color = c("green3", "grey70"),
      breaks = c(-0.5, 0.5, 1.5),
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      treeheight_row = 0,
      annotation_col = anno,
      fontsize_row = 8,
      fontsize_col = 8,
      legend_breaks = c(0, 1),
      legend_labels = c("Present", "Missing"),
      annotation_colors = mycolors
    )
  } else {
    missingnessHeatmap <- FALSE
  }

  ## Potential Specificity ----
  # Add groups back to miss.df
  miss.df$Group = TheData$Group

  # For columns which have any amount missing, find those with p.val < 0.05
  anyMissing <- colnames(miss.df[, which(missing < 0.9)])

  # Create df with only significant metabolites
  # res.spec <- NULL
  # for(i in 1:length(anyMissing)){
  #   metabolite <- anyMissing[i]
  #
  #   #formula of metabolite ~ Group and then run through glm
  #   f <- as.formula(paste0(metabolite,"~ Group"))
  #   fit <- glm(f,data = miss.df, method = "brglmFit")
  #
  #   res.t <- data.frame(Compound = metabolite, summary(emmeans(fit,pairwise~Group)$contrasts))
  #   res.spec <- rbind(res.spec, res.t[res.t$p.value < 0.05,])
  # }

  ## Remove Missing ----
  # Table showing what will be removed
  missingnessTable <- names(missing[missing < 0.75])

  # Remove the missing
  TheData <- TheData[, !colnames(TheData) %in% missingnessTable]

  return(list(
    TheData = TheData,
    missingnessHeatmap = missingnessHeatmap,
    missingnessTable = missingnessTable
  ))
}

#' QRILC imputation
#'
#' @param dataSet.mvs input data?
#' @param tune.sigma Tuning param passed to tmvtnorm::rtmvnorm sigma
#' @returns A list of results
#'
#' @export
impute.QRILC <- function(dataSet.mvs, tune.sigma = 1) {
  nFeatures = dim(dataSet.mvs)[1]
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  QR.obj = list()
  for (i in 1:nSamples) {
    curr.sample = dataSet.mvs[, i]
    pNAs = length(which(is.na(curr.sample))) / length(curr.sample)
    upper.q = 0.99
    q.normal = qnorm(
      seq(
        (pNAs + 0.001),
        (upper.q + 0.001),
        (upper.q - pNAs) / (upper.q * 100)
      ),
      mean = 0,
      sd = 1
    )
    q.curr.sample = quantile(
      curr.sample,
      probs = seq(0.001, (upper.q + 0.001), 0.01),
      na.rm = T
    )
    temp.QR = lm(q.curr.sample ~ q.normal)
    QR.obj[[i]] = temp.QR
    mean.CDD = temp.QR$coefficients[1]
    sd.CDD = as.numeric(temp.QR$coefficients[2])
    data.to.imp = tmvtnorm::rtmvnorm(
      n = nFeatures,
      mean = mean.CDD,
      sigma = sd.CDD * tune.sigma,
      upper = qnorm(
        (pNAs +
          0.001),
        mean = mean.CDD,
        sd = sd.CDD
      ),
      algorithm = c("gibbs")
    )
    curr.sample.imputed = curr.sample
    curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(
      curr.sample
    ))]
    dataSet.imputed[, i] = curr.sample.imputed
  }
  results = list(dataSet.imputed, QR.obj)
  return(results)
}
