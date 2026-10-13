# Get metabolites where proportion present is less than a set threshold

Get metabolites where proportion present is less than a set threshold

## Usage

``` r
IdentifyMissing(qcData, presence_threshold = 0.75)
```

## Arguments

- qcData:

  QC'd dataset

- presence_threshold:

  A threshold for non-missing data below which a metabolite should be
  dropped

## Value

Vector of metabolite names below the threshold
