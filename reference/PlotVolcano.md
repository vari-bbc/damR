# Plot volcano for a single contrast

Plot volcano for a single contrast

## Usage

``` r
PlotVolcano(
  resultTable,
  contrastGroups,
  contrastTitle = NULL,
  alphaLevel = 0.05,
  cap = NULL,
  ...
)
```

## Arguments

- resultTable:

  The `topTable` results from
  [`pAdjustData()`](https://vari-bbc.github.io/damR/reference/pAdjustData.md)

- contrastGroups:

  Two contrast group names

- contrastTitle:

  Optional title for the Volcano plot

- alphaLevel:

  Significance cutoff level

- cap:

  Absolute log fold-change cap to cap the `logFC` x axis so that
  `-cap < logFC < cap`.

- ...:

  Geom aesthetics passed to
  [`VolcanoPlotter()`](https://vari-bbc.github.io/damR/reference/VolcanoPlotter.md).

## Value

Volcano ggplot
