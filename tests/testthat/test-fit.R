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
  exdata <- system.file(package = "metaboapp") |>
    file.path("extdata") |>
    list.files("example_data.csv", full.names = TRUE) |>
    ReadFile()
  qc_data <- QCData(exdata)$TheData
  complete_data <- RemoveMissingness(qc_data, IdentifyMissing(qc_data, 0.75))
  complete_data[, 7:ncol(complete_data)] <- impute.QRILC(
    complete_data[, 7:ncol(complete_data)]
  )[[1]]
  expect_no_error(
    ebayes.fit <- FitTheData(
      complete_data,
      c("AMW20p_AMW20q", "AMW20p_AMW20s", "BDp_AMW20s"),
      "title"
    )
  )
  expect_s4_class(ebayes.fit, "MArrayLM")
})
