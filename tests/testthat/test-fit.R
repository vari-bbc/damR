test_that("MakeTheContrasts", {
  expect_equal(
    MakeTheContrasts(LETTERS[1:4], "contrast_title"),
    c(
      contrast_title1 = "A - B",
      contrast_title2 = "A - C",
      contrast_title3 = "A - D",
      contrast_title4 = "B - C",
      contrast_title5 = "B - D",
      contrast_title6 = "C - D"
    )
  )
  expect_equal(
    MakeTheContrasts(LETTERS[1:4], "contrast_title", "Group"),
    c(
      contrast_title1 = "GroupA - GroupB",
      contrast_title2 = "GroupA - GroupC",
      contrast_title3 = "GroupA - GroupD",
      contrast_title4 = "GroupB - GroupC",
      contrast_title5 = "GroupB - GroupD",
      contrast_title6 = "GroupC - GroupD"
    )
  )
})

test_that("FitTheData", {
  exdata <- ReadFile(exdata_path)
  qc_data <- QCData(exdata)
  complete_data <- RemoveMissing(qc_data)
  imputed_data <- impute.QRILC(complete_data)
  expect_no_error(
    ebayes.fit <- DiffMetabolites(complete_data)
  )
})
