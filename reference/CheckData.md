# Check that data meets requirements

Checks that the values are numeric and positive, and that `Sample`
column is present and none of its identifiers have more than 1 dash
character

## Usage

``` r
CheckData(input_data)
```

## Arguments

- input_data:

  `data.frame` of the input data from `read.csv`.

## Value

input_data `data.frame` of the input data if no errors are found
