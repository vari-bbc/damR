test_that("ReadFile", {
  read_in_data <- ReadFile(exdata_path)
  expected_colnames <- c(
    "Sample",
    "SampleID",
    "Group",
    "Subgroup1",
    "Subgroup2",
    "Batch",
    "Metabolite",
    "Type",
    "value"
  )
  expect_equal(colnames(read_in_data), expected_colnames)

  # test that #-starting lines are ignored
  comment_file <- file.path(extdata, "comment_data.csv")
  comment_filtered <- ReadFile(comment_file)
  expect_all_false(grepl("^#", comment_filtered$Sample))
})

test_that("CheckData", {
  # extra dash in one sample id

  tidy <- input_data[, 1:2]
  colnames(tidy) <- c("Sample", "value")
  more_dashes <- tidy
  more_dashes[1, 1] <- sub("_", "-", more_dashes[1, 1])

  # character in data column
  not_numeric <- tidy
  not_numeric[1, 2] <- "a"

  expect_error(
    CheckData(more_dashes),
    "Each 'Sample' column entry may contain only one dash \\(max of two conditions\\)"
  )
  expect_error(
    CheckData(not_numeric),
    "The data columns must all be numeric"
  )
})
