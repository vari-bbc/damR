exdata <- ReadFile(exdata_path)
qc_data <- exdata
# qc_data <- QCData(exdata)$TheData

n <- 100
dummy.df <- data.table("Metabolite" = rep(c("HighMissing", "LowMissing"), 1, each = n))
dummy.df$value <- seq_len(nrow(dummy.df))

dummy.df <- rbind(
  # 0.5 missing
  dummy.df[Metabolite == "HighMissing"][seq(ceiling(n * 0.50)), value := NA],
  # 0.26 missing
  dummy.df[Metabolite == "LowMissing"][seq(ceiling(n * 0.25)), value := NA]
)

test_that("IdentifyMissing", {
  expect_equal(IdentifyMissing(dummy.df, 0.80), c("HighMissing", "LowMissing"))
  missing_cols <- IdentifyMissing(dummy.df, 0.75)
  expect_equal(missing_cols, c("HighMissing"))
  removed <- RemoveMissing(dummy.df, 0.75)
  expect_equal(removed[, unique("LowMissing")], "LowMissing")
})
