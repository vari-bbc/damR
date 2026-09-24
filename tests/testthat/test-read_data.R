test_that("ReadFile", {
  read_in_data <- ReadFile(exdata_path)
  expected_colnames <- c("Sample", "SampleID", "Group", "Subgroup1", "Subgroup2", "Batch")
  expect_equal(colnames(read_in_data)[1:6], expected_colnames)

  # test that #-starting lines are ignored
  comment_file <- file.path(extdata, "comment_data.csv")
  comment_unfiltered <- read.csv(comment_file)
  comment_filtered <- ReadFile(comment_file)
  expect_equal(nrow(comment_unfiltered), 5) # test of the test file to ensure row count
  expect_equal(nrow(comment_filtered), 4)
  expect_all_false(grepl("^#", comment_filtered$Sample))
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
