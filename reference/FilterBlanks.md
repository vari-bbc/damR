# Blank filter data by batch

This function will set filtered values to NA. To get a table of the
blank values that get dropped use
[`GetBlankDropTable()`](https://vari-bbc.github.io/damR/reference/GetBlankDropTable.md)
instead, and then `dropBlankFiltered` to do the filtering.

## Usage

``` r
FilterBlanks(inputData, blankFoldThreshold = 5, blankPattern = "PB|blank")
```

## Arguments

- inputData:

  Input metabolite data from
  [`ReadFile()`](https://vari-bbc.github.io/damR/reference/ReadFile.md)

- blankFoldThreshold:

  Fold change threshold below which data are considered blanks

- blankPattern:

  Sample name pattern to use to capture blank. Will match
  case-insensitively on substring using
  [`data.table::like`](https://rdrr.io/pkg/data.table/man/like.html)
  (`%ilike`)

## Value

inputData without the metabolites that don't meet the threshold and
without the blanks used for filtering

## Details

Adapted from the proteomics blank_filter_by_batch() in "Updated Blank
Filter.R". For each metabolite, a sample value must exceed the loudest
process blank in its own batch by at least `blankFoldThreshold` fold;
values that don't are set to NA. The blanks are only used to decide what
to drop – surviving values keep their original intensities, they are not
blank-subtracted. Expects TheData metabolite columns to already be log2
transformed (as they are by the time QCData runs), which makes the fold
change a simple subtraction.
