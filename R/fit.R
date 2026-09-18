#' Produce contrasts from comparisons
#'
#' @param comparisonSelections Vector of groups to compare
#' @param comparisonsTitle Name prefix for the contrasts entries, suffixed by their index
#' @param prefix Optional prefix to add to each contrast level. Used for limma contrasts
#' @returns A named vector with contrasts in the form A - B
#'
#' @importFrom utils combn
#'
#' @export
MakeTheContrasts <- function(comparisonSelections, comparisonsTitle, prefix = "") {
  res <- combn(paste0(prefix, comparisonSelections), 2) |>
    apply(2, paste, collapse = " - ")
  names(res) <- paste0(comparisonsTitle, seq(length(res)))
  res
}

#' Fit model with limma
#'
#' @param cleanData The QC'd and imputed `data.frame`
#' @param theContrasts Vector of groups names to contrast
#' @param contrastTitle A title for the contrasts groups used for plot titles
#' @param contrastLevel The column from the data containing contrasts levels
#' @param repeated_measures Were there repeated measures?
#' @returns ebayes fit from limma
#'
#' @importFrom limma contrasts.fit duplicateCorrelation eBayes lmFit topTable
#' @importFrom stats model.matrix p.adjust
#'
#' @export
FitTheData <- function(
  cleanData,
  theContrasts,
  contrastTitle,
  contrastLevel = "Group",
  repeated_measures = "no"
) {
  colnames(cleanData) <- make.names(colnames(cleanData))
  meta_cols <- c("Sample", "SampleID", "Group", paste0("Subgroup", 1:2), "Batch")
  metaboliteCols <- which(colnames(cleanData) %notin% meta_cols)

  cleanData$contrastLevel <- cleanData[, contrastLevel]

  design <- model.matrix(~ 0 + contrastLevel, data = cleanData) # ERW Make model
  num_data <- data.frame(t(cleanData[, metaboliteCols])) # ERW Remove non metabolite cols and flip dataframe
  print(num_data[1:5, 1:8])

  # Get named vector with contrasts named by the title given
  contrasts.design <- MakeTheContrasts(theContrasts, contrastTitle, prefix = "contrastLevel")

  # lmFit the metabolite data using the model above
  set.seed(111)
  if (repeated_measures == "yes") {
    # if repeated, then you need a duplicate correlation in your lmfit
    block <- cleanData$Replicate
    duplCorl <- duplicateCorrelation(num_data, design, block = block)

    fit <- lmFit(
      num_data,
      design,
      block = block,
      correlation = duplCorl$consensus.correlation
    )
  } else {
    fit <- lmFit(num_data, design)
  }

  # Make the contrasts to use with our fit using our contrasts and the model
  cont <- limma::makeContrasts(contrasts = contrasts.design, levels = design)
  fit2 <- contrasts.fit(fit, cont) # Get estimated coefficients and standard errors
  eBayes(fit2, trend = TRUE, robust = TRUE) # Get t-stats, f-stat, and log odds
}


#' Mulitple testing correction
#'
#' @param ebayesFit The fit from FitTheData
#' @param theContrasts The contrasts for comparison
#' @returns P-adjusted fit
#'
pAdjustData <- function(ebayesFit, theContrasts) {
  # Create succinct table to p-adjust
  adjustedData <- NULL
  for (i in 1:length(unlist(theContrasts))) {
    compoundTable <- topTable(
      ebayesFit,
      adjust.method = "none",
      p.value = 1,
      number = Inf,
      coef = i
    )
    compounds <- rownames(compoundTable)
    adjustedData <- rbind(
      adjustedData,
      data.frame(
        Compound = compounds,
        Contrast = unlist(theContrasts)[i],
        subset(compoundTable, select = c("logFC", "adj.P.Val"))
      ),
      row.names = NULL
    )
  }

  # Do the p-adjustment
  adjustedData$adj.P.Val <- p.adjust(adjustedData$adj.P.Val, method = "BH")

  return(adjustedData)
}
