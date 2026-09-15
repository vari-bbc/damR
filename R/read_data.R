#' Read and process input data file
#'
#' @param filePath Path to the input file
#' @returns Processed data.frame
#'
#' @importFrom dplyr relocate
#'
#' @export
ReadFile <- function(filePath) {
  theData <- read.csv(filePath, check.names = FALSE)

  ## Pull sample info ----
  # Set up for seperating sample into info
  theData$SampleID <- NA
  theData$Group <- NA
  theData$Subgroup1 <- NA
  theData$Subgroup2 <- NA
  theData$Batch <- NA
  theData = theData |>
    relocate(SampleID, Group, Subgroup1, Subgroup2, Batch, .after = 'Sample')

  # Fill out info
  i <- 1
  for (i in 1:length(theData$Sample)) {
    # Seperate sample by splitting by "_"
    theSample <- theData$Sample[i]
    theSplit <- unlist(strsplit(theSample, "_"))

    # First is sampleID, last is batch
    theData[i, "SampleID"] <- theSplit[1]
    theData[i, "Batch"] <- theSplit[length(theSplit)]

    # Combine middle and split by "-" to get subgroup
    thegroups <- theSplit[2:(length(theSplit) - 1)]
    allGroups <- paste0(thegroups, collapse = "_")
    groupSplit <- unlist(strsplit(allGroups, "-"))

    # Group is allGroups with the dash as an underscore
    allGroups <- gsub("-", "_", allGroups)
    theData[i, "Group"] <- allGroups

    # Subgroups are from the group split with a max of 2
    theData[i, "Subgroup1"] <- groupSplit[1]
    theData[i, "Subgroup2"] <- groupSplit[2]
  }

  return(theData)
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
