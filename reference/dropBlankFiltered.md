# Drop data that fails to meet the blank fold change threshold as found by `GetBlankDropTable()`

Performs a join on the table of value to be dropped and sets the value
column of the match to NA, returning only non-blank data

## Usage

``` r
dropBlankFiltered(inputData, dropTable, blanks)
```

## Arguments

- inputData:

  Input metabolite data from
  [`ReadFile()`](https://vari-bbc.github.io/metaboapp/reference/ReadFile.md)

- dropTable:

  `data.table` of data points that don't meet the blank fold threshold
  from
  [`GetBlankDropTable()`](https://vari-bbc.github.io/metaboapp/reference/GetBlankDropTable.md)

- blanks:

  Vector of Sample names of blank samples

## Value

filtered `inputData` where values in `dropTable` are set to NA and blank
samples removed
