# Plot PCA

PCA scatter plot with blanks and QC pool observations marked with a
border

## Usage

``` r
PlotPCA(
  pca,
  x = "PC1",
  y = "PC2",
  colorBy = "Group",
  pointSize = 3,
  pointAlpha = 0.8
)
```

## Arguments

- pca:

  The PCA and variance explained data from
  [`GetPCA()`](https://vari-bbc.github.io/damR/reference/GetPCA.md)

- x, y:

  The principal components (prefixed with "PC") to plot on the X and Y
  axes

- colorBy:

  The column by which to color the points (using `aes(fill = ...)`)

- pointSize:

  The size of the
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- pointAlpha:

  The transparency of the
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

## Value

A named list where `pca` is the PCA data column bound to `fullData` and
`varex` which contains the percent variance explained for each PC.
