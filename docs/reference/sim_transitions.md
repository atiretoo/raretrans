# Simulate population projection matrices

`sim_transition` generates a list of simulated population projection
matrices from the provided parameters and prior distributions.

## Usage

``` r
sim_transitions(
  TF,
  N,
  P = NULL,
  alpha = 1e-05,
  beta = 1e-05,
  priorweight = -1,
  samples = 1
)
```

## Arguments

- TF:

  A list of two matrices, T and F, as ouput by
  [`projection.matrix`](https://rdrr.io/pkg/popbio/man/projection.matrix.html).

- N:

  A vector of observed stages at start of transition.

- P:

  A matrix of the priors for each column. Defaults to uniform.

- alpha:

  A matrix of the prior parameter for each stage. Impossible stage
  combinations marked with NA_real\_.

- beta:

  A matrix of the prior parameter for each stage. Impossible stage
  combinations marked with NA_real\_.

- priorweight:

  total weight for each column of prior as a percentage of sample size
  or 1 if negative

- samples:

  The number of matrices to return.

## Value

Always returns a list.

## Examples

``` r
# Build a simple 3-stage TF list
T_mat <- matrix(c(0.5, 0.3, 0.0,
                  0.2, 0.4, 0.1,
                  0.0, 0.1, 0.7), nrow = 3, ncol = 3)
F_mat <- matrix(c(0.0, 0.0, 1.5,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0), nrow = 3, ncol = 3)
TF <- list(T = T_mat, F = F_mat)
N  <- c(10, 5, 8)

# Simulate 10 population projection matrices using default uninformative priors
set.seed(42)
mats <- sim_transitions(TF, N, samples = 10)
#> Warning: length(alpha | beta) != order^2: only using first value of alpha and beta
length(mats)   # 10 matrices
#> [1] 10
mats[[1]]      # first simulated matrix
#>           [,1]         [,2]       [,3]
#> [1,] 0.7601157 0.0098559821 0.00239957
#> [2,] 0.1753338 0.5190143831 0.13319806
#> [3,] 1.4237391 0.0001148738 0.76055177
```
