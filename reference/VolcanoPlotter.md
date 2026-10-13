# Volcano plot

Volcano plot

## Usage

``` r
VolcanoPlotter(
  pd,
  alphaLevel,
  pointColors = c(`Up regulated` = "#276419", `Down regulated` = "#8E0152", NS = "grey"),
  pointSize = 2,
  alphaLineColor = "red",
  alphaLineType = "dotted",
  alphaLineWidth = 0.5,
  interactive = FALSE
)
```

## Arguments

- pd:

  Plot data produced by
  [`VolcanoData()`](https://vari-bbc.github.io/metaboapp/reference/VolcanoData.md)

- alphaLevel:

  Significance cutoff level

- pointColors:

  Named vector of colors for `"Up regulated"`, `"Down regulated"`, and
  `"NS"` points

- pointSize:

  Size of the
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

- alphaLineColor, alphaLineType, alphaLineWidth:

  Color, type and width of the line used to mark the significance
  cutoff. See `?ggplot2::geom_hline()`.

- interactive:

  Whether to plot as interactive with
  [`ggiraph::geom_point_interactive()`](https://davidgohel.github.io/ggiraph/reference/geom_point_interactive.html)

## Value

Volcano ggplot
