exdata_path <- system.file(package = "metaboapp") |>
  file.path("extdata") |>
  list.files("example_data.csv", full.names = TRUE)
input_data <- read.csv(exdata_path)
