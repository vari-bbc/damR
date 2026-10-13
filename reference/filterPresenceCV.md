# Identify the highest non-NA proportion with lowest coefficient of variance

Identify the highest non-NA proportion with lowest coefficient of
variance

## Usage

``` r
filterPresenceCV(pPresent, cv, naTolerance)
```

## Arguments

- pPresent:

  Non-NA proportions for a metabolite across methods

- cv:

  Coefficient of variance for a metabolite across methods

- naTolerance:

  The absolute acceptable non-NA proportion difference below with to
  pick the lowest `cv`.

## Value

The index of the row with least NAs or, if non-NA proportions are within
`naTolerance`, the row with lowest CV

## Details

This function requires that the input `pPresent` vector is sorted in
descending order. The between proportion differences are classified as
being above or below the `naTolerance`. If the first difference, always
between the top two highest proportions when inputs are sorted, is
larger than `naTolerance`, it returns 1 since the biggest `pPresent`
exceeds the tolerance. Otherwise it returns the lowest coefficient of
variance among the values with difference less than `naTolerance`. The
indexing with [`which()`](https://rdrr.io/r/base/which.html) works since
`diff` returns a vector with one less than the input.
