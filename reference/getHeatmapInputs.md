# Create the data and annotations for heatmaps

Create the data and annotations for heatmaps

## Usage

``` r
getHeatmapInputs(dt, anno.pal)
```

## Arguments

- dt:

  The long-form metabolite data from qc for missingness or imputed data
  that's subset to significant metabolites.

- anno.pal:

  RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation
  colors

## Value

list with Data matrix, annotation data frame, and colors for the
subgroups

## Details

Used for plotting missingness and the significant metabolite heatmap.
