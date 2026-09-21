#' Read and process input data file
#'
#' @details
#' Splits the Sample column by `"_"` into SampleID, Group and Batch columns.
#' Then splits the Group column Subgroup1 and Subgroup2 on `"-"`. The `"-"` in
#' the Group column is replaced with `"_"` to prevent issues with making
#' contrasts that need a `"-"` to delineate between contrast groups.
#'
#' @param filePath Input file path
#' @returns Processed data.frame
#'
#' @importFrom rlang .data
#'
#' @export
ReadFile <- function(filePath) {
  read.csv(filePath) |>
    CheckData() |>
    tidyr::separate("Sample", into = c("SampleID", "Group", "Batch"), sep = "_", remove = FALSE) |>
    tidyr::separate("Group", into = c("Subgroup1", "Subgroup2"), sep = "-", remove = FALSE) |>
    dplyr::mutate("Group" = gsub("-", "_", .data$Group))
}

#' Check that the data contains the required columns
#'
#' Checks that the Sample column is present and none of its identifiers have
#' more than 1 dash character
#'
#' @param input_data `data.frame` of the input data from `read.csv`.
#' @returns input_data `data.frame` of the input data if no errors are found
#'
#' @export
#' @keywords internal
CheckData <- function(input_data) {
  stopifnot("The data must contain a 'Sample' column." = "Sample" %in% colnames(input_data))
  stopifnot(
    "Each 'Sample' column entry may contain only one dash (max of two conditions)" = all(
      stringr::str_count(input_data$Sample, "-") == 1
    )
  )
  input_data
}
