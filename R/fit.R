#' Produce contrasts from comparisons
#'
#' @param labels Vector of groups to compare
#' @param title Vector base name for the contrasts entries. Will be suffixed by index
#' @param prefix Optional prefix to add to each contrast level. Used for limma contrasts
#' @returns A named vector with contrasts in the form A - B
#'
#' @importFrom stats as.formula
#' @importFrom utils combn
#'
#' @export
#' @examples
#' # untitled
#' MakeTheContrasts(letters[1:4])
#' # with names
#' MakeTheContrasts(letters[1:4], "name")
#' # with names and prefix; used for model matrix creation
#' MakeTheContrasts(letters[1:4], "name", "Group")
MakeTheContrasts <- function(labels, title = "", prefix = "") {
  res <- combn(paste0(prefix, labels), 2) |>
    apply(2, paste, collapse = " - ")
  names(res) <- paste0(title, seq(length(res)))
  res
}

#' Fit model with limma
#'
#' @param cleanData The QC'd and imputed `data.frame`
#' @param theContrasts Vector of groups names to contrast. Defaults to the
#' unique elements of the `Group` column.
#' @param contrastTitle A title for the contrasts groups used for plot titles.
#' Defaults to `""`.
#' @param contrastLevel The column from the data containing contrasts levels
#' @param repeated_measures Were there repeated measures?
#' @returns ebayes fit from limma
#'
#' @importFrom collapse funique num_vars
#' @importFrom limma contrasts.fit duplicateCorrelation eBayes lmFit topTable
#' @importFrom stats model.matrix p.adjust
#'
#' @export
FitTheData <- function(
  cleanData,
  theContrasts = NULL,
  contrastTitle = NULL,
  contrastLevel = "Group",
  repeated_measures = "no"
) {
  theContrasts <- theContrasts %||% funique(cleanData[[contrastLevel]])
  contrastTitle <- contrastTitle %||% ""

  mat <- fselect(cleanData, c("Sample", "Group", "Metabolite", "value")) |>
    pivot(ids = c("Sample", "Group"), names = "Metabolite", how = "wider")
  colnames(cleanData) <- make.names(colnames(cleanData))

  design <- model.matrix(as.formula(paste("~ 0 +", contrastLevel)), data = mat) # ERW Make model
  num_data <- t(num_vars(mat))
  print(num_data[1:5, 1:8])

  # Get named vector with contrasts named by the title given
  contrasts.design <- MakeTheContrasts(theContrasts, contrastTitle, prefix = contrastLevel)

  # lmFit the metabolite data using the model above
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


#' Multiple testing correction within contrasts
#'
#' @param ebayesFit The fit from FitTheData
#' @param theContrasts Vector of groups to compare
#' @param contrastTitle Name prefix for the contrasts entries, suffixed by their index
#'
#' @returns `data.frame` of topTable hits with adjusted P-values and logFC
#'
#' @importFrom collapse fselect
#'
#' @export
pAdjustData <- function(ebayesFit, theContrasts = NULL, contrastTitle = NULL) {
  contrastLabels <- MakeTheContrasts(theContrasts, contrastTitle)

  tt_list <- lapply(seq_along(contrastLabels), function(i) {
    tt <- topTable(ebayesFit, adjust.method = "BH", p.value = 1, number = Inf, coef = i)
    tt$Contrast <- contrastLabels[i]
    tt$Compound <- rownames(tt)
    rownames(tt) <- NULL
    tt
  })
  do.call(rbind, tt_list) |>
    fselect(c("Compound", "Contrast", "logFC", "adj.P.Val"))
}
