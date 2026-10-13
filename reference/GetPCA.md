# Compute PCA and variance explained

Compute PCA and variance explained

## Usage

``` r
GetPCA(fullData)
```

## Arguments

- fullData:

  The full dataset including non-numeric data

## Value

A named list where `pca` is the PCA data column bound to `fullData` and
`varex` which contains the percent variance explained for each PC.
