test_that("MakeTheContrasts", {
  expect_equal(
    MakeTheContrasts(LETTERS[1:4], "contrast_title"),
    list(contrast_title = c("A - B", "A - C", "A - D", "B - C", "B - D", "C - D"))
  )
  expect_equal(
    MakeTheContrasts(LETTERS[1:4], "contrast_title", "Group"),
    list(
      contrast_title = c(
        "GroupA - GroupB",
        "GroupA - GroupC",
        "GroupA - GroupD",
        "GroupB - GroupC",
        "GroupB - GroupD",
        "GroupC - GroupD"
      )
    )
  )
})

test_that("FitTheData", {
  exdata <- system.file(package = "metabbclomics") |>
    file.path("extdata") |>
    list.files("example_data.csv", full.names = TRUE) |>
    ReadFile()
  qc_data <- QCData(exdata)$TheData
  complete_data <- RemoveMissingness(qc_data)[[1]]
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
