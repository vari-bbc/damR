# Read and process input data file

Read and process input data file

## Usage

``` r
ReadFile(filePath)
```

## Arguments

- filePath:

  Input file path

## Value

Processed data.frame

## Details

Any lines starting with '#' are ignored to avoid reading template
formatting. Pivots the data from wide form into a tidy/long format.
Splits the Sample column by `"_"` into SampleID, Group and Batch
columns. Then splits the Group column Subgroup1 and Subgroup2 on `"-"`.
The `"-"` in the Group column is replaced with `"_"` to prevent issues
with making contrasts that need a `"-"` to delineate between contrast
groups. The Metabolite column is split on `"_"` to identify the
extraction method. All values \<= 0 are set to NA and other values log2
transformed. The following requirements must be met for this function to
work correctly:

- Sample column data in the form "SampleID_Subgroup1-Subgroup2_Batch"

- Metabolite column data in the form "MetaboliteName_ExtractionMethod"
