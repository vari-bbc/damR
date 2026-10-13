# Multiple testing correction within contrasts. `DiffMetabolites()` should be preferred to using this function directly.

Multiple testing correction within contrasts.
[`DiffMetabolites()`](https://vari-bbc.github.io/damR/reference/DiffMetabolites.md)
should be preferred to using this function directly.

## Usage

``` r
pAdjustData(ebayesFit, theContrasts, contrastLevel = "Group")
```

## Arguments

- ebayesFit:

  The fit from FitTheData

- theContrasts:

  Vector of contrasts specified in the format "GroupA - GroupB". The
  "Group" string, specified by `contrastLevel` is stripped.

- contrastLevel:

  The column from the data containing contrasts levels. `model.matrix`
  adds this column name as a prefix to each group so this is passed in
  with the prefix to
  [`FitTheData()`](https://vari-bbc.github.io/damR/reference/FitTheData.md)
  but can now be stripped contrast names in the result.

## Value

`data.frame` of topTable hits with adjusted P-values and logFC
