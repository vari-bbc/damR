# Produce contrasts from comparisons

Produce contrasts from comparisons

## Usage

``` r
MakeTheContrasts(labels, title = "", prefix = "")
```

## Arguments

- labels:

  Vector of groups to compare

- title:

  Vector base name for the contrasts entries. Will be suffixed by index

- prefix:

  Optional prefix to add to each contrast level. Used for limma
  contrasts

## Value

A named vector with contrasts in the form A - B

## Examples

``` r
# untitled
MakeTheContrasts(letters[1:4])
#>       1       2       3       4       5       6 
#> "a - b" "a - c" "a - d" "b - c" "b - d" "c - d" 
# with names
MakeTheContrasts(letters[1:4], "name")
#>   name1   name2   name3   name4   name5   name6 
#> "a - b" "a - c" "a - d" "b - c" "b - d" "c - d" 
# with names and prefix; used for model matrix creation
MakeTheContrasts(letters[1:4], "name", "Group")
#>             name1             name2             name3             name4 
#> "Groupa - Groupb" "Groupa - Groupc" "Groupa - Groupd" "Groupb - Groupc" 
#>             name5             name6 
#> "Groupb - Groupd" "Groupc - Groupd" 
```
