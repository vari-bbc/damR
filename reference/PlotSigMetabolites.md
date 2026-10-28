# Plot heatmap of significantly differentially abundant metabolites

Plot heatmap of significantly differentially abundant metabolites

## Usage

``` r
PlotSigMetabolites(
  diffMetabolites,
  imputedData,
  fcThreshold = 2,
  alpha = 0.05,
  scale = FALSE,
  anno.pal = c("Set1", "Set2"),
  ...
)
```

## Arguments

- diffMetabolites:

  The differentially abundant metabolite data from
  [`DiffMetabolites()`](https://vari-bbc.github.io/damR/reference/DiffMetabolites.md)
  used to select metabolites for the heatmap

- imputedData:

  Imputed metabolite data from
  [`impute.QRILC()`](https://vari-bbc.github.io/damR/reference/impute.QRILC.md)
  used to produce the heatmap

- fcThreshold:

  The fold change to filter metabolites as 'significant'

- alpha:

  Significance level to filter the differential abundance p-value

- scale:

  Whether to scale the imputed data matrix to mean 0 and standard
  deviation of 1

- anno.pal:

  RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation
  colors

- ...:

  Arguments passed on to
  [`ComplexHeatmap::pheatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/pheatmap.html)

  `color`

  :   The same as in `pheatmap`. Here you don't necessarily need to
      generate a long color vector. The discrete colors sent to
      [`colorRampPalette`](https://rdrr.io/r/grDevices/colorRamp.html)
      are also OK here. E.g.
      `colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)`
      can be simply replaced as
      `rev(brewer.pal(n = 7, name = "RdYlBu"))`.

  `kmeans_k`

  :   The same as in `pheatmap`.

  `breaks`

  :   The same as in `pheatmap`.

  `border_color`

  :   The same as in `pheatmap`.

  `cellwidth`

  :   The same as in `pheatmap`.

  `cellheight`

  :   The same as in `pheatmap`.

  `cluster_rows`

  :   The same as in `pheatmap`.

  `cluster_cols`

  :   The same as in `pheatmap`.

  `clustering_distance_rows`

  :   The same as in `pheatmap`.

  `clustering_distance_cols`

  :   The same as in `pheatmap`.

  `clustering_method`

  :   The same as in `pheatmap`.

  `clustering_callback`

  :   The same as in `pheatmap`.

  `cutree_rows`

  :   The same as in `pheatmap`.

  `cutree_cols`

  :   The same as in `pheatmap`.

  `treeheight_row`

  :   The same as in `pheatmap`.

  `treeheight_col`

  :   The same as in `pheatmap`.

  `legend`

  :   The same as in `pheatmap`.

  `legend_breaks`

  :   The same as in `pheatmap`.

  `legend_labels`

  :   The same as in `pheatmap`.

  `annotation_row`

  :   The same as in `pheatmap`.

  `annotation`

  :   The same as in `pheatmap`.

  `annotation_legend`

  :   The same as in `pheatmap`.

  `annotation_names_row`

  :   The same as in `pheatmap`.

  `annotation_names_col`

  :   The same as in `pheatmap`.

  `drop_levels`

  :   Enforced to be `TRUE`.

  `show_rownames`

  :   The same as in `pheatmap`.

  `show_colnames`

  :   The same as in `pheatmap`.

  `main`

  :   The same as in `pheatmap`.

  `fontsize`

  :   The same as in `pheatmap`.

  `fontsize_row`

  :   The same as in `pheatmap`.

  `fontsize_col`

  :   The same as in `pheatmap`.

  `angle_col`

  :   The same as in `pheatmap`.

  `display_numbers`

  :   The same as in `pheatmap`.

  `number_format`

  :   The same as in `pheatmap`.

  `number_color`

  :   The same as in `pheatmap`.

  `fontsize_number`

  :   The same as in `pheatmap`.

  `gaps_row`

  :   The same as in `pheatmap`.

  `gaps_col`

  :   The same as in `pheatmap`.

  `labels_row`

  :   The same as in `pheatmap`.

  `labels_col`

  :   The same as in `pheatmap`.

  `filename`

  :   Not supported.

  `width`

  :   Not supported.

  `height`

  :   Not supported.

  `silent`

  :   Not supported.

  `na_col`

  :   The same as in `pheatmap`.

  `name`

  :   Name of the heatmap. This argument is passed to
      [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

  `fontfamily`

  :   Font family for row and column names.

  `fontfamily_row`

  :   Font family for row names.

  `fontfamily_col`

  :   Font family for column names.

  `fontface`

  :   Font face for row and column names.

  `fontface_row`

  :   Font face for row names.

  `fontface_col`

  :   Font face for column names.

  `heatmap_legend_param`

  :   Pass to
      [`Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

  `run_draw`

  :   Whether to run `draw()` function to the heatmap object.

## Value

ComplexHeatmap::pheatmap

## Details

Additional plot aesthetics can be controlled through arguments passed to
[`ComplexHeatmap::pheatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/pheatmap.html)
though `...`. However, `mat`, `annotation_col`, and `annotation_colors`
are set by `damR` so modifying those arguments may fail. Use `anno.pal`
to modify the annotation colors.
