# Fit model with limma across all. `DiffMetabolites()` should be preferred to using this function directly.

Fit model with limma across all.
[`DiffMetabolites()`](https://vari-bbc.github.io/damR/reference/DiffMetabolites.md)
should be preferred to using this function directly.

## Usage

``` r
FitTheData(
  cleanData,
  theContrasts,
  contrastLevel = "Group",
  repeated_measures = "no"
)
```

## Arguments

- cleanData:

  The QC'd and imputed `data.frame`

- theContrasts:

  Vector of contrasts specified in the format "GroupA - GroupB".
  Produced by
  [`MakeTheContrasts()`](https://vari-bbc.github.io/damR/reference/MakeTheContrasts.md)
  with the `prefix` argument set to "Group" by default as thats the
  column usually used for contrasts.

- contrastLevel:

  The column from the data containing `groups` levels, usually
  `"Group"`.

- repeated_measures:

  Were there repeated measures?

## Value

ebayes fit from limma
