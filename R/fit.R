#' Produce contrasts from comparisons
#'
#' @param comparisonSelections Vector of groups to compare
#' @param comparisonsTitle Name for the contrasts
#' @returns A named list with contrasts in the form A-B
#'
#' @importFrom utils combn
#'
#' @export
MakeTheContrasts <- function(comparisonSelections, comparisonsTitle) {
  res <- combn(comparisonSelections, 2) |>
    apply(2, paste, collapse = "-") |>
    list(... = _)
  names(res) <- comparisonsTitle
  res
}


#' Fit model with limma
#'
#' @param TheData The cleaned? data
#' @param theContrasts The contrasts for comparison
#' @returns ebayes fit from limma
#'
#' @importFrom limma contrasts.fit duplicateCorrelation eBayes lmFit topTable
#' @importFrom stats model.matrix p.adjust
#'
#' @export
FitTheData <- function(TheData, theContrasts) {
  metaboliteCols <- 7:ncol(TheData)

  colnames(TheData) <- make.names(colnames(TheData))

  TheData$contrastLevel <- TheData[, contrastLevel]

  # Make model
  # ERW
  design <- model.matrix(~ 0 + contrastLevel, data = TheData)

  # Remove non metabolite cols and flip dataframe
  # ERW
  ys <- data.frame(t(TheData[, metaboliteCols]))

  print(ys[1:5, 1:8])

  # Get named vector with contrasts named by the title given
  contrasts.design <- sapply(strsplit(unlist(theContrasts), " "), function(x) {
    x[nchar(x) > 1] = paste0("contrastLevel", x[nchar(x) > 1])
    return(paste(x, collapse = " "))
  })

  # lmFit the metabolite data using the model above
  set.seed(111)
  if (repeated_measures == "yes") {
    # if repeated, then you need a duplicate correlation in your lmfit
    block <- TheData$Replicate
    duplCorl <- duplicateCorrelation(ys, design, block = block)

    fit <- lmFit(
      ys,
      design,
      block = block,
      correlation = duplCorl$consensus.correlation
    )
  } else {
    fit <- lmFit(ys, design)
  }

  # Make the contrasts to use with our fit using our contrasts and the model
  cont <- limma::makeContrasts(contrasts = contrasts.design, levels = design)

  # Get estimated coefficients and standard errors
  set.seed(111)
  fit2 <- contrasts.fit(fit, cont)

  # Get t-stats, f-stat, and log odds
  set.seed(111)
  ebayesFit <- eBayes(fit2, trend = TRUE, robust = TRUE)

  return(ebayesFit)
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
