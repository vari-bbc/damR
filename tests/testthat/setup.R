library("data.table")
extdata <- system.file(package = "damR") |> file.path("extdata")
exdata_path <- list.files(extdata, "example_data.csv", full.names = TRUE)
input_data <- read.csv(exdata_path)
