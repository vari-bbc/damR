#' Blank filter data by batch
#'
#' This function will set filtered values to NA. To get a table of the blank
#' values that get dropped use `getDropTable()` instead, and then
#' `dropBlankFiltered` to do the filtering.
#'
#' @details
#' Adapted from the proteomics blank_filter_by_batch() in "Updated Blank
#' Filter.R". For each metabolite, a sample value must exceed the loudest
#' process blank in its own batch by at least `blankFoldThreshold` fold; values
#' that don't are set to NA. The blanks are only used to decide what to drop --
#' surviving values keep their original intensities, they are not
#' blank-subtracted. Expects TheData metabolite columns to already be log2
#' transformed (as they are by the time QCData runs), which makes the fold
#' change a simple subtraction.
#'
#' @param inputData Input metabolite data from `ReadFile()`
#' @param blankFoldThreshold Fold change threshold below which data are
#' considered blanks
#' @param blankPattern Sample name pattern to use to capture blank. Will match case-insensitively on substring using `data.table::like` (`%ilike`)
#'
#' @returns inputData without the metabolites that don't meet the threshold and
#' without the blanks used for filtering
#'
#' @importFrom collapse %=% .c
#'
#' @export
BlankFilter <- function(inputData, blankFoldThreshold = 5, blankPattern = "PB|blank") {
  # get list with values to drop and the blanks
  dropTable <- blanks <- NULL
  .c(dropTable, blanks) %=% getDropTable(inputData, blankFoldThreshold, blankPattern)
  if (is.null(dropTable)) {
    return(inputData)
  }
  dropBlankFiltered(inputData, dropTable, blanks)
}

#' Drop data that fails to meet the blank fold change threshold as found by
#' `getDropTable()`
#'
#' Performs a join on the table of value to be dropped and sets the value
#' column of the match to NA, returning only non-blank data
#'
#' @param inputData Input metabolite data from `ReadFile()`
#' @param dropTable `data.table` of data points that don't meet the blank fold threshold from `getDropTable()`
#' @param blanks Vector of Sample names of blank samples
#'
#' @returns filtered `inputData` where values in `dropTable` are set to NA and
#' blank samples removed
#'
#' @importFrom data.table copy
#'
#' @export
dropBlankFiltered <- function(inputData, dropTable, blanks) {
  Sample <- value <- NULL
  copy(inputData)[dropTable, on = c("Sample", "Metabolite", "Batch"), value := NA] |>
    fsubset(Sample %notin% blanks)
}

#' Identify data that falls under the provided threshold by batch across batches
#'
#' @inherit BlankFilter details
#'
#' @inheritParams BlankFilter
#'
#' @returns named list where `dropTable` is a data.table with Sample,
#' Metabolite, value, and Batch and `blanks` are the Sample names corresponding
#' to blank samples. This table does not contain any of the blanks or other
#' control sequences.
#'
#' @importFrom collapse fsubset funique
#' @importFrom data.table %like% %ilike%
#'
#' @export
getDropTable <- function(inputData, blankFoldThreshold = 5, blankPattern = "PB|blank") {
  Sample <- NULL
  blanks <- funique(fsubset(inputData, Sample %ilike% "blank")[["Sample"]])
  if (length(blanks) == 0) {
    message(sprintf("No blanks were found in this dataset using the pattern: '%s'", blankPattern))
    return(list(dropTable = NULL, blanks = NULL))
  }
  list(
    # split by Batch and blank filter on each batch data.table, referenced by .SD
    dropTable = inputData[, getBatchDropTable(.SD, log2(blankFoldThreshold), blanks), by = "Batch"],
    blanks = blanks
  )
}

#' Identify data that falls under the threshold for one batch
#'
#' Helper for BlankFilter and getDropTable
#'
#' @inherit BlankFilter details
#'
#' @param batchData A single batch of input data. This must contain a blank
#' @param logFoldThreshold log2 of fold change threshold below which data are
#' considered blanks
#' @param blanks Vector of Sample names of blank samples
#'
#' @returns data.table with Sample, Metabolite and value where the value does
#' not meet the blank filtering threshold
#'
#' @importFrom collapse fgroup_by fmutate fsubset fsubset fsummarize join
#' @importFrom data.table data.table .SD
#'
#' @export
#' @keywords internal
getBatchDropTable <- function(batchData, logFoldThreshold, blanks) {
  value <- Sample <- Metabolite <- max_blank <- diff <- NULL

  non_na_blanks <- fsubset(batchData, Sample %in% blanks & !is.na(value)) |>
    fgroup_by(Metabolite) |>
    fsummarize(max_blank = max(value))

  if (nrow(non_na_blanks) == 0) {
    return(data.table())
  }

  fsubset(batchData, Sample %notin% blanks) |>
    join(non_na_blanks, on = "Metabolite", verbose = 0) |>
    fmutate(diff = value - max_blank) |>
    fsubset(diff < logFoldThreshold)
}
