# Get differential metabolites with log fold-changes and adjusted pvalues based on input comparison groups

Fits model on the full input dataset with limma

## Usage

``` r
DiffMetabolites(
  cleanData,
  groups = NULL,
  contrastLevel = "Group",
  repeated_measures = "no"
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

## Value

LogFC and adjusted P-values of
[`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html) for
each contrast
