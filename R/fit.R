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

#' Get differential metabolites with log fold-changes and adjusted pvalues
#' based on input comparison groups
#'
#' Fits model on the full input dataset with limma
#'
#' @param cleanData The QC'd and imputed `data.frame`
#' @param groups Vector of groups names to contrast. Defaults to the
#' unique elements of the `Group` column but the column used can be specified
#' using `contrastLevel`.
#' @param contrastLevel The column from the data containing `groups` levels,
#' usually `"Group"`.
#' @param repeated_measures Were there repeated measures?
#'
#' @returns LogFC and adjusted P-values of `limma::topTable()` for each contrast
#'
#' @export
DiffMetabolites <- function(
  cleanData,
  groups = NULL,
  contrastLevel = "Group",
  repeated_measures = "no"
) {
  groups <- groups %||% funique(cleanData[[contrastLevel]])
  theContrasts <- MakeTheContrasts(groups, prefix = contrastLevel)
  ebayesFit <- FitTheData(
    cleanData,
    theContrasts,
    contrastLevel,
    repeated_measures
  )
  pAdjustData(ebayesFit, theContrasts, contrastLevel)
}

#' Fit model with limma across all. `DiffMetabolites()` should be preferred to
#' using this function directly.
#'
#' @inheritParams DiffMetabolites cleanData contrastLevel repeated_measures
#' @param theContrasts Vector of contrasts specified in the format "GroupA -
#' GroupB". Produced by `MakeTheContrasts()` with the `prefix` argument set to
#' "Group" by default as thats the column usually used for contrasts.
#'
#' @returns ebayes fit from limma
#'
#' @importFrom collapse funique num_vars
#' @importFrom limma contrasts.fit duplicateCorrelation eBayes lmFit makeContrasts topTable
#' @importFrom stats model.matrix p.adjust
#'
#' @export
FitTheData <- function(
  cleanData,
  theContrasts,
  contrastLevel = "Group",
  repeated_measures = "no"
) {
  mat <- fselect(cleanData, c("Sample", "Group", "Metabolite", "value")) |>
    pivot(ids = c("Sample", "Group"), names = "Metabolite", how = "wider")
  colnames(cleanData) <- make.names(colnames(cleanData))

  design <- model.matrix(as.formula(paste("~ 0 +", contrastLevel)), data = mat) # ERW Make model
  num_data <- t(num_vars(mat))

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
  cont <- makeContrasts(contrasts = theContrasts, levels = design)
  fit2 <- contrasts.fit(fit, cont) # Get estimated coefficients and standard errors
  eBayes(fit2, trend = TRUE, robust = TRUE) # Get t-stats, f-stat, and log odds
}


#' Multiple testing correction within contrasts. `DiffMetabolites()` should be
#' preferred to using this function directly.
#'
#' @param ebayesFit The fit from FitTheData
#' @param theContrasts Vector of contrasts specified in the format "GroupA -
#' GroupB". The "Group" string, specified by `contrastLevel` is stripped.
#' @param contrastLevel The column from the data containing contrasts levels.
#' `model.matrix` adds this column name as a prefix to each group so this is
#' passed in with the prefix to `FitTheData()` but can now be stripped contrast
#' names in the result.
#'
#' @returns `data.frame` of topTable hits with adjusted P-values and logFC
#'
#' @importFrom collapse fselect
#'
#' @export
pAdjustData <- function(ebayesFit, theContrasts, contrastLevel = "Group") {
  # extract just the contrasts of interest from the full fit
  contrastLabels <- stringi::stri_replace_all_fixed(theContrasts, contrastLevel, "")

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
