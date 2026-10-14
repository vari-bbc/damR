# Select the point geom based on interactivity

Select the point geom based on interactivity

## Usage

``` r
getGeomPoint(interactive)
```

## Arguments

- interactive:

  Whether the point should be interactive

## Value

[`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
if interactive else
[`ggiraph::geom_point_interactive()`](https://davidgohel.github.io/ggiraph/reference/geom_point_interactive.html)
