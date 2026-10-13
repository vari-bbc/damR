# Run QC on input data

Run QC on input data

## Usage

``` r
QCData(
  inputData,
  filterBlanks = TRUE,
  filterExtraction = TRUE,
  blankFoldThreshold = 5,
  naTolerance = 0.05,
  blankPattern = "^PB|blank",
  controlPattern = "^PB|^TC|^QC"
)
```

## Arguments

- inputData:

  Input metabolite data

- filterBlanks:

  Whether to run blank filtering. See
  [`FilterBlanks()`](https://vari-bbc.github.io/damR/reference/FilterBlanks.md).

- filterExtraction:

  Whether to filter for the best extraction method. See
  [`FilterExtraction()`](https://vari-bbc.github.io/damR/reference/FilterExtraction.md).

- blankFoldThreshold:

  The blank fold change threshold. The difference between every
  metabolite and its batch's blank must be at least this value or higher
  to remain in the dataset.

- naTolerance:

  The absolute acceptable non-NA proportion difference for a given
  metabolite across separation methods like IP, HILICA and T3Acidic.
  Metabolites with proportion differences under this value will be
  picked based on their coefficient of variance. See
  [`FilterExtraction()`](https://vari-bbc.github.io/damR/reference/FilterExtraction.md).

- blankPattern:

  String pattern to denote blank samples

- controlPattern:

  String pattern to denote control samples like QC

## Details

For each metabolite's separation method, it selects the one with the
lowest proportion NAs, or if similar proportions, the lowest coefficient
of variance. See
[`FilterExtraction()`](https://vari-bbc.github.io/damR/reference/FilterExtraction.md).
Then removes

- data points below a fold change threshold compared to blank samples.
  See `BlankFilter()`.

- metabolites with zero variance

- Blank samples (after blank filtering)

- Control samples
