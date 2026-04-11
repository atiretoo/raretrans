# Helper functions for generating different priors

Extract the number of individuals in each stage from a dataframe of
transitions.

## Usage

``` r
get_state_vector(transitions, stage = NULL, sort = NULL)
```

## Arguments

- transitions:

  a dataframe of observations of individuals in different stages

- stage:

  the name of the variable with the stage information

- sort:

  a vector of stage names in the desired order. Default is the order of
  levels in stage.

## Value

a vector of the counts of observations in each level of stage.

## Examples

``` r
data("L_elto")

# Extract one population at one census period
onepop <- L_elto[L_elto$POPNUM == 250 & L_elto$year == 5, ]
onepop$stage <- factor(onepop$stage, levels = c("p", "j", "a"))

# Count individuals per stage in the order p, j, a
get_state_vector(onepop, stage = "stage", sort = c("p", "j", "a"))
#> [1] 11 47 34
```
