#' Read and process input data file
#'
#' @details
#' Splits the Sample column by `"_"` into SampleID, Group and Batch columns.
#' Then splits the Group column Subgroup1 and Subgroup2 on `"-"`. The `"-"` in
#' the Group column is replaced with `"_"` to prevent issues with making
#' contrasts that need a `"-"` to delineate between contrast groups.
#'
#' @param filePath Input file path @returns Processed data.frame
#'
#' @importFrom rlang .data
#'
#' @export
ReadFile <- function(filePath) {
  read.csv(filePath) |>
    tidyr::separate("Sample", into = c("SampleID", "Group", "Batch"), sep = "_", remove = FALSE) |>
    tidyr::separate("Group", into = c("Subgroup1", "Subgroup2"), sep = "-", remove = FALSE) |>
    dplyr::mutate("Group" = gsub("-", "_", .data$Group))
}

#' Check that the data contains the required columns
#'
#' @param theData Data read by ReadFile
#' @returns error message to be passed onto Shiny; FALSE if no errors
#'
#' @export
#' @keywords internal
CheckData <- function(theData) {
  errorText <- "Data Errors:\n"
  isErrors <- FALSE

  # Check if the data has "Sample" column
  if (!"Sample" %in% colnames(theData)) {
    errorText <- paste0(errorText, "The data must contain a 'Sample' column.\n")
    isErrors <- TRUE
  }

  # Check if the sample column has more than one dash in each entry
  if (any(stringr::str_count(theData$Sample, "-") > 1)) {
    errorText <- paste0(
      errorText,
      "Each entry in the 'Sample' column must have a max of one dash (max of two conditions).\n"
    )
    isErrors <- TRUE
  }

  #if any errors, return the error text, otherwise return the data
  if (isErrors) {
    return(errorText)
  } else {
    return(FALSE)
  }
}
