# Select the best extraction method for each metabolite

For each metabolite's separation method, it selects the one with the
lowest proportion NAs, or the lowest coefficient of variance when the
proportions are similar.

## Usage

``` r
FilterExtraction(inputData, naTolerance)
```

## Arguments

- inputData:

  Input metabolite data

- naTolerance:

  The absolute acceptable non-NA proportion difference for a given
  metabolite across separation methods like IP, HILICA and T3Acidic. The
  extraction method choice for metabolites with proportion differences
  under this value is based on their coefficient of variance. At 0.05,
  non-NA proportions of 0.95 and 0.93 would be considered equivalent and
  the metabolite would be selected based on the lowest coefficient of
  variance instead.

## Value

`data.table`, The `inputData` without metabolites from extraction
methods with high NA-proportion and variation
