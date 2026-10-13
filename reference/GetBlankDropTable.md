# Identify data that falls under the provided threshold by batch across batches

Identify data that falls under the provided threshold by batch across
batches

## Usage

``` r
GetBlankDropTable(inputData, blankFoldThreshold = 5, blankPattern = "PB|blank")
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

named list where `dropTable` is a data.table with Sample, Metabolite,
value, and Batch and `blanks` are the Sample names corresponding to
blank samples. This table does not contain any of the blanks or other
control sequences.

## Details

Adapted from the proteomics blank_filter_by_batch() in "Updated Blank
Filter.R". For each metabolite, a sample value must exceed the loudest
process blank in its own batch by at least `blankFoldThreshold` fold;
values that don't are set to NA. The blanks are only used to decide what
to drop – surviving values keep their original intensities, they are not
blank-subtracted. Expects TheData metabolite columns to already be log2
transformed (as they are by the time QCData runs), which makes the fold
change a simple subtraction.
