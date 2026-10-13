# Identify data that falls under the threshold for one batch

Helper for FilterBlanks and GetBlankDropTable

## Usage

``` r
getBatchDropTable(batchData, logFoldThreshold, blanks)
```

## Arguments

- batchData:

  A single batch of input data. This must contain a blank

- logFoldThreshold:

  log2 of fold change threshold below which data are considered blanks

- blanks:

  Vector of Sample names of blank samples

## Value

data.table with Sample, Metabolite and value where the value does not
meet the blank filtering threshold

## Details

Adapted from the proteomics blank_filter_by_batch() in "Updated Blank
Filter.R". For each metabolite, a sample value must exceed the loudest
process blank in its own batch by at least `blankFoldThreshold` fold;
values that don't are set to NA. The blanks are only used to decide what
to drop – surviving values keep their original intensities, they are not
blank-subtracted. Expects TheData metabolite columns to already be log2
transformed (as they are by the time QCData runs), which makes the fold
change a simple subtraction.
