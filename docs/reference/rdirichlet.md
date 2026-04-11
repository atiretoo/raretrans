# Dirichlet distributed random numbers

Generate n random vectors distributed according to a Dirichlet
distribution. Each row of the returned matrix is a random vector that
sums to 1.

## Usage

``` r
rdirichlet(n, alpha)
```

## Source

copied from package `MCMCpack` to avoid a dependency. That code was
taken from Greg's Miscellaneous Functions (gregmisc). His code was based
on code posted by Ben Bolker to R-News on 15 Dec 2000.

## Arguments

- n:

  The number of random vectors to generate

- alpha:

  A vector of parameters

## Value

The function returns a matrix with n rows and `length(alpha)` columns

## Examples

``` r
# Generate 5 random probability vectors from a Dirichlet(1, 2, 3) distribution
set.seed(42)
rdirichlet(5, c(1, 2, 3))
#>             [,1]      [,2]      [,3]
#> [1,] 0.358883528 0.1643524 0.4767640
#> [2,] 0.008081945 0.4629716 0.5289465
#> [3,] 0.540750715 0.2614354 0.1978138
#> [4,] 0.276963518 0.3291203 0.3939162
#> [5,] 0.343456431 0.1620774 0.4944662
```
