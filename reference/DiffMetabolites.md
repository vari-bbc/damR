# Get differential metabolites with log fold-changes and adjusted pvalues based on input comparison groups

Fits model on the full input dataset with limma

## Usage

``` r
DiffMetabolites(
  cleanData,
  groups = NULL,
  contrastLevel = "Group",
  repeated_measures = "no",
  reverseContrasts = TRUE
)
```

## Arguments

- cleanData:

  The QC'd and imputed `data.frame`

- groups:

  Vector of groups names to contrast. Defaults to the unique elements of
  the `Group` column but the column used can be specified using
  `contrastLevel`.

- contrastLevel:

  The column from the data containing `groups` levels, usually
  `"Group"`.

- repeated_measures:

  Were there repeated measures?

- reverseContrasts:

  Make contrasts in both directions. This makes it easy to examine both
  A - B and B - A. If set to false, the contrast order depends on the
  order of the input groups and the combinations produced by
  [`MakeTheContrasts()`](https://vari-bbc.github.io/damR/reference/MakeTheContrasts.md)
  which may not be available in the output when reversed.

## Value

LogFC and adjusted P-values of
[`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html) for
each contrast
