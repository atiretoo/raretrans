# Combine prior and data for transition matrix

`fill_transition` returns the expected value of the transition matrix
combining observed transitions for one time step and a prior

## Usage

``` r
fill_transitions(TF, N, P = NULL, priorweight = -1, returnType = "T")
```

## Arguments

- TF:

  A list of two matrices, T and F, as ouput by
  [`projection.matrix`](https://rdrr.io/pkg/popbio/man/projection.matrix.html).

- N:

  A vector of observed transitions.

- P:

  A matrix of the priors for each column. Defaults to uniform.

- priorweight:

  total weight for each column of prior as a percentage of sample size
  or 1 if negative

- returnType:

  A character vector describing the desired return value. Defaults to
  "T" the transition matrix.

## Value

The return value depends on parameter returnType.

- A - the summed matrix "filled in" using a dirichlet prior

- T - just the filled in transition matrix

- TN - the augmented matrix of fates – use in calculating the CI or for
  simulation

## Examples

``` r
# Build a simple 3-stage TF list (transition + fertility matrices)
T_mat <- matrix(c(0.5, 0.3, 0.0,
                  0.2, 0.4, 0.1,
                  0.0, 0.1, 0.7), nrow = 3, ncol = 3)
F_mat <- matrix(c(0.0, 0.0, 1.5,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0), nrow = 3, ncol = 3)
TF <- list(T = T_mat, F = F_mat)
N  <- c(10, 5, 8)

# Default: return filled transition matrix T
fill_transitions(TF, N)
#>            [,1]      [,2]       [,3]
#> [1,] 0.47727273 0.2083333 0.02777778
#> [2,] 0.29545455 0.3750000 0.11666667
#> [3,] 0.02272727 0.1250000 0.65000000

# Return the full population matrix A = T + F
fill_transitions(TF, N, returnType = "A")
#>           [,1]      [,2]       [,3]
#> [1,] 0.4772727 0.2083333 0.02777778
#> [2,] 0.2954545 0.3750000 0.11666667
#> [3,] 1.5227273 0.1250000 0.65000000

# Use a prior weight equal to the sample size
fill_transitions(TF, N, priorweight = 1)
#>       [,1]  [,2]  [,3]
#> [1,] 0.375 0.225 0.125
#> [2,] 0.275 0.325 0.175
#> [3,] 0.125 0.175 0.475
```
