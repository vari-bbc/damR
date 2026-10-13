# Helper to set up annotation colors for Sample groups in heatmap

Helper to set up annotation colors for Sample groups in heatmap

## Usage

``` r
.getAnnoColors(anno.df, anno.pal)
```

## Arguments

- anno.df:

  Annotation `data.frame` containing Subgroups as columns

- anno.pal:

  RColorBrewer palette to use for Subgroup1 and Subgroup2 annotation
  colors

## Value

Named list of named vectors of colors. List names are `anno.df` column
names and vector names are of the unique elements in `anno.df`.

## Examples

``` r
(anno.df <- data.frame(Subgroup1 = rep(letters[1:2], 3), Subgroup2 = letters[3:8]))
#>   Subgroup1 Subgroup2
#> 1         a         c
#> 2         b         d
#> 3         a         e
#> 4         b         f
#> 5         a         g
#> 6         b         h
.getAnnoColors(anno.df, c("Set1", "Set2"))
#> Warning: minimal value for n is 3, returning requested palette with 3 different levels
#> $Subgroup1
#>         a         b 
#> "#E41A1C" "#377EB8" 
#> 
#> $Subgroup2
#>         c         d         e         f         g         h 
#> "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" 
#> 
```
