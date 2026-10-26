#' Run QC on input data
#'
#' @details
#' For each metabolite's separation method, it selects the one with the
#' lowest proportion NAs, or if similar proportions, the lowest coefficient
#' of variance. See `FilterExtraction()`.
#' Then removes
#'   - data points below a fold change threshold compared to blank samples. See
#'   `BlankFilter()`.
#'   - metabolites with zero variance
#'   - Blank samples (after blank filtering)
#'   - Control samples
#'
#' @param inputData Input metabolite data
#' @param filterExtraction Whether to filter for the best extraction method.
#' See `FilterExtraction()`.
#' @param filterBlanks Whether to run blank filtering. See `FilterBlanks()`.
#' @param naTolerance The absolute acceptable non-NA proportion difference for a
#' given metabolite across separation methods like IP, HILICA and T3Acidic.
#' Metabolites with proportion differences under this value will be picked
#' based on their coefficient of variance. See `FilterExtraction()`.
#' @param blankFoldThreshold The blank fold change threshold. The difference
#' between every metabolite and its batch's blank must be at least this value
#' or higher to remain in the dataset.
#' @param blankPattern String pattern to denote blank samples
#' @param controlPattern String pattern to denote control samples like QC
#'
#' @returns QC'd `data.table`
#'
#' @importFrom collapse fgroup_by fsubset fsummarise
#' @importFrom stats sd
#'
#' @export
QCData <- function(
  inputData,
  filterBlanks = TRUE,
  filterExtraction = TRUE,
  blankFoldThreshold = 5,
  naTolerance = 0.05,
  blankPattern = "^PB|blank",
  controlPattern = "^PB|^TC|^QC"
) {
  if (filterExtraction) {
    inputData <- FilterExtraction(inputData, naTolerance)
  }

  if (filterBlanks) {
    inputData <- FilterBlanks(inputData, blankFoldThreshold, blankPattern)
  }

  zeroVarMetabolites <- getZeroVar(inputData)
  zeroVarDropped <- dropMetabolites(zeroVarMetabolites, inputData)

  Sample <- NULL
  fsubset(zeroVarDropped, !(Sample %like% controlPattern))
}

#' Identify metabolites with zero variance, including any that are all NA
#'
#' @param inputData Input metabolite data
#'
#' @returns Vector of metabolites with zero variance
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
    fsubset(zv | is.na(zv)) |>
    _[["Metabolite"]]
}
