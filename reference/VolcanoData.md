# Produce volcano plot data for a single contrast

Produce volcano plot data for a single contrast

## Usage

``` r
VolcanoData(
  resultTable,
  contrastGroups,
  contrastTitle,
  alphaLevel = 0.05,
  cap = 10
)
```

## Arguments

- resultTable:

  The `topTable` results from
  [`pAdjustData()`](https://vari-bbc.github.io/metaboapp/reference/pAdjustData.md)

- contrastGroups:

  Two contrast group names

- contrastTitle:

  Optional title for the Volcano plot

- alphaLevel:

  Significance cutoff level

- cap:

  Absolute log fold-change cap to cap the `logFC` x axis so that
  `-cap < logFC < cap`.

## Value

data.table ready for plotting with
[`VolcanoPlotter()`](https://vari-bbc.github.io/metaboapp/reference/VolcanoPlotter.md)
