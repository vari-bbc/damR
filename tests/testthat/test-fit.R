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
