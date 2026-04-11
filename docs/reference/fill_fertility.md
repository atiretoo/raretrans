# Combine prior and data for fertility matrix

`fill_fertility` returns the expected value of the fertility matrix
combining observed recruits for one time step and a Gamma prior for each
column.

## Usage

``` r
fill_fertility(
  TF,
  N,
  alpha = 1e-05,
  beta = 1e-05,
  priorweight = -1,
  returnType = "F"
)
```

## Arguments

- TF:

  A list of two matrices, T and F, as ouput by
  [`projection.matrix`](https://rdrr.io/pkg/popbio/man/projection.matrix.html).

- N:

  A vector of observed stage distribution.

- alpha:

  A matrix of the prior parameter for each stage. Impossible stage
  combinations marked with NA_real\_.

- beta:

  A matrix of the prior parameter for each stage. Impossible stage
  combinations marked with NA_real\_.

- priorweight:

  total weight for each column of prior as a percentage of sample size
  or 1 if negative

- returnType:

  A character vector describing the desired return value. Defaults to
  "F" the fertility matrix

## Value

The return value depends on parameter returnType.

- A - the summed matrix "filled in" using a Gamma prior

- F - just the filled in fertility matrix

- ab - the posterior parameters alpha and beta as a list.

## Details

Assumes that only one stage reproduces ... needs generalizing.

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

# Only adults (stage 3) reproduce; mark non-reproducing entries with NA
alpha_mat <- matrix(c(NA, NA, 0.5,
                      NA, NA, NA,
                      NA, NA, NA), nrow = 3, ncol = 3)
beta_mat  <- matrix(c(NA, NA, 1.0,
                      NA, NA, NA,
                      NA, NA, NA), nrow = 3, ncol = 3)

# Default: return filled fertility matrix F
fill_fertility(TF, N, alpha = alpha_mat, beta = beta_mat)
#>          [,1] [,2] [,3]
#> [1,] 0.000000    0    0
#> [2,] 0.000000    0    0
#> [3,] 1.409091    0    0

# Return the full population matrix A = T + F
fill_fertility(TF, N, alpha = alpha_mat, beta = beta_mat, returnType = "A")
#>          [,1] [,2] [,3]
#> [1,] 0.500000  0.2  0.0
#> [2,] 0.300000  0.4  0.1
#> [3,] 1.409091  0.1  0.7

# Return the posterior alpha and beta parameters
fill_fertility(TF, N, alpha = alpha_mat, beta = beta_mat, returnType = "ab")
#> $alpha
#>      [,1] [,2] [,3]
#> [1,]   NA   NA   NA
#> [2,]   NA   NA   NA
#> [3,] 15.5   NA   NA
#> 
#> $beta
#>      [,1] [,2] [,3]
#> [1,]   NA   NA   NA
#> [2,]   NA   NA   NA
#> [3,]   11   NA   NA
#> 
```
