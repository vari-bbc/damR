# Plot missingness heatmap

Selects the necessary columns to make a wide form data frame. Makes the
presence/absence matrix with the Sample column and numeric columns and
an annotation data frame for annotation colors using the Subgroups.

## Usage

``` r
PlotMissingness(qcData, anno.pal = c("Set1", "Set2"), showRownames = TRUE)
```

## Arguments

- qcData:

  Input data that has been QC'd

- anno.pal:

  RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation
  colors

- showRownames:

  Whether to show metabolite names as rownames in the heatmap

## Value

ComplexHeatmap::pheatmap
