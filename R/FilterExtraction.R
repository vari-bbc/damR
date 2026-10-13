#' Select the best extraction method for each metabolite
#'
#' For each metabolite's separation method, it selects the one with the
#' lowest proportion NAs, or the lowest coefficient of variance when the
#' proportions are similar.
#'
#' @param inputData Input metabolite data
#' @param naTolerance The absolute acceptable non-NA proportion difference for a
#' given metabolite across separation methods like IP, HILICA and T3Acidic.
#' The extraction method choice for metabolites with proportion differences
#' under this value is based on their coefficient of variance. At 0.05, non-NA
#' proportions of 0.95 and 0.93 would be considered equivalent and the
#' metabolite would be selected based on the lowest coefficient of variance
#' instead.
#'
#' @returns `data.table`, The `inputData` without metabolites from extraction
#' methods with high NA-proportion and variation
#'
#' @importFrom collapse fgroup_by fmean fsummarize
#' @importFrom data.table data.table .SD setorder
#'
#' @export
FilterExtraction <- function(inputData, naTolerance) {
  value <- pPresent <- cv <- Metabolite <- NULL
  prop.present <- fgroup_by(inputData, c("Metabolite", "Type")) |>
    fsummarize(
      pPresent = 1 - fmean(is.na(value)),
      cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE)
    ) |>
    setorder(Metabolite, -pPresent)

  # run filterPresenceCV on each Metabolite group, represented by .SD
  selection <- prop.present[,
    .SD[filterPresenceCV(pPresent, cv, naTolerance), "Type"],
    by = "Metabolite"
  ]
  join(inputData, selection, on = c("Metabolite", "Type"), how = "inner", verbose = 0)
}

#' Identify the highest non-NA proportion with lowest coefficient of variance
#'
#' @details
#' This function requires that the input `pPresent` vector is sorted in
#' descending order. The between proportion differences are classified as being
#' above or below the `naTolerance`. If the first difference, always between
#' the top two highest proportions when inputs are sorted, is larger than
#' `naTolerance`, it returns 1 since the biggest `pPresent` exceeds the
#' tolerance. Otherwise it returns the lowest coefficient of variance among the
#' values with difference less than `naTolerance`. The indexing with `which()`
#' works since `diff` returns a vector with one less than the input.
#'
#' @param pPresent Non-NA proportions for a metabolite across methods
#' @param cv Coefficient of variance for a metabolite across methods
#' @param naTolerance The absolute acceptable non-NA proportion difference
#' below with to pick the lowest `cv`.
#'
#' @returns The index of the row with least NAs or, if non-NA proportions are
#' within `naTolerance`, the row with lowest CV
#'
#' @export
#' @keywords internal
filterPresenceCV <- function(pPresent, cv, naTolerance) {
  isLargeDiff <- abs(diff(pPresent)) > naTolerance
  firstLargeDiff <- which(isLargeDiff)[1]
  ifelse(
    is.na(firstLargeDiff), # if all are FALSE or within naTolerance
    which.min(cv),
    which.min(cv[seq(firstLargeDiff)]) # always 1 if 1st if larger than naTolerance
  )
}
