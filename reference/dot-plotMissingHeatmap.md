# Helper to plot missingness heatmap

Helper to plot missingness heatmap

## Usage

``` r
.plotMissingHeatmap(bin.mat, anno.df, anno.colors, showRownames)
```

## Arguments

- bin.mat:

  Binary matrix where 0 is missing and 1 is present

- anno.df:

  Annotation data.frame with the sample groups

- anno.colors:

  List named with `anno.df` column names containing named color vectors
  where the color vector names correspond to the elements of anno. See
  [`.getAnnoColors()`](https://vari-bbc.github.io/metaboapp/reference/dot-getAnnoColors.md)

- showRownames:

  Whether to show metabolite names as rownames in the heatmap

## Value

ComplexHeatmap
