#' Read and process input data file
#'
#' @details
#' Splits the Sample column by `"_"` into SampleID, Group and Batch columns.
#' Then splits the Group column Subgroup1 and Subgroup2 on `"-"`. The `"-"` in
#' the Group column is replaced with `"_"` to prevent issues with making
#' contrasts that need a `"-"` to delineate between contrast groups. Any lines
#' starting with '#' are ignored to avoid reading template formatting. All
#' values <= 0 are set to NA and other values log2 transformed.
#'
#' @param filePath Input file path
#' @returns Processed data.frame
#'
#' @importFrom data.table fcase fread %like% melt setcolorder tstrsplit :=
#'
#' @export
ReadFile <- function(filePath) {
  x <- value <- NULL
  df.wide <- fread(filePath, comment.char = "#", integer64 = "double")
  df.tidy <- suppressWarnings(melt(df.wide, id.vars = "Sample", variable.name = "Metabolite"))
  df.tidy[value == 0, value := NA] # drop 0 so no log2(0)

  CheckData(df.tidy)

  # split up Sample column; `env` specifies columns to prevent R CMD check NOTES
  df.tidy[, c("SampleID", "Group", "Batch") := tstrsplit(x, "_"), env = list(x = "Sample")] |>
    _[, paste0("Subgroup", 1:2) := tstrsplit(x, "-"), env = list(x = "Group")] |>
    _[, "Group" := gsub("-", "_", x), env = list(x = "Group")] |>
    _[,
      # get IP vs HILICA
      "type" := fcase(
        x %like% "_IP$"     , "IP"     ,
        x %like% "_HILICA$" , "HILICA" ,
        default = NA
      ),
      env = list(x = "Metabolite")
    ]
  # return log2 values
  df.tidy[, "value" := log2(x), env = list(x = "value")]
  cols <- c("Sample", "SampleID", "Group", "Subgroup1", "Subgroup2", "Batch", "Metabolite", "type", "value")
  setcolorder(df.tidy, cols)
  df.tidy[]
}

#' Check that data meets requirements
#'
#' Checks that the values are numeric and positive, and that `Sample` column is
#' present and none of its identifiers have more than 1 dash character
#'
#' @param input_data `data.frame` of the input data from `read.csv`.
#' @returns input_data `data.frame` of the input data if no errors are found
#'
#' @export
#' @keywords internal
CheckData <- function(input_data) {
  stopifnot(
    "Each 'Sample' column entry may contain only one dash (max of two conditions)" = all(
      stringi::stri_count_fixed(input_data$Sample, "-") == 1
    )
  )
  stopifnot("The data columns must all be numeric" = is.numeric(input_data$value))
  stopifnot("The data columns must all be positive" = all(input_data$value >= 0, na.rm = TRUE))
  invisible(NULL)
}
