exdata <- ReadFile(exdata_path)[, 1:10]
qc_data <- QCData(exdata)$TheData
nSamples <- nrow(qc_data)
dummy_data <- seq(nSamples)
dummy.df <- data.frame(A = dummy_data, B = dummy_data, C = dummy_data)
qc_data <- cbind(qc_data[, 1:6], dummy.df)

# 0.5 missing
qc_data[, "A"][seq(ceiling(nSamples * 0.50))] <- NA
# 0.26 missing
qc_data[, "B"][seq(ceiling(nSamples * 0.26))] <- NA

test_that("IdentifyMissing", {
  missing_cols <- IdentifyMissing(qc_data, 0.75)
  expect_equal(missing_cols, c("A", "B"))
  removed <- RemoveMissingness(qc_data, missing_cols)
  expect_true(all(colnames(removed) %notin% c("A", "B")))
})
