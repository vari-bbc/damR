test_that("ReadFile", {
  input_data <- ReadFile(exdata_path)
  expected_colnames <- c("Sample", "SampleID", "Group", "Subgroup1", "Subgroup2", "Batch")
  expect_equal(colnames(input_data)[1:6], expected_colnames)
})

test_that("CheckData", {
  # no sample column
  missing_sample <- input_data[, -1]

  # extra dash in one sample id
  more_dashes <- input_data
  more_dashes[1, 1] <- sub("_", "-", more_dashes[1, 1])

  # character in data column
  not_numeric <- input_data
  not_numeric[1, 2] <- "a"

  expect_error(CheckData(missing_sample), "The data must contain a 'Sample' column.")
  expect_error(
    CheckData(more_dashes),
    "Each 'Sample' column entry may contain only one dash \\(max of two conditions\\)"
  )
  expect_error(
    CheckData(not_numeric),
    "The data columns must all be numeric"
  )
})
