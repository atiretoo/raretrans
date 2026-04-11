# Calculate beta credible intervals for all transition matrix entries

Computes the marginal posterior beta credible intervals for every entry
in the transition matrix, including the probability of dying. The
marginal posterior distribution of each transition probability follows a
beta distribution derived from the Dirichlet-multinomial model, using
the augmented fate matrix (TN) returned by
[`fill_transitions`](https://atiretoo.github.io/raretrans/reference/fill_transitions.md).

## Usage

``` r
transition_CrI(
  TF,
  N,
  P = NULL,
  priorweight = -1,
  ci = 0.95,
  stage_names = NULL
)
```

## Arguments

- TF:

  A list of two matrices, T and F, as output by
  [`projection.matrix`](https://rdrr.io/pkg/popbio/man/projection.matrix.html).

- N:

  A vector of observed stage distribution at the start of the transition
  period.

- P:

  A matrix of priors for each column. Defaults to uniform.

- priorweight:

  Total weight for each column of prior as a percentage of sample size,
  or 1 if negative. Defaults to -1 (uninformative).

- ci:

  Credible interval width as a probability between 0 and 1. Defaults to
  0.95 (95% credible interval).

- stage_names:

  Optional character vector of stage names in the same order as the
  columns of the transition matrix. If `NULL`, names are taken from
  `colnames(TF$T)`, or generic labels `"Stage 1"`, `"Stage 2"`, etc. are
  used.

## Value

A data frame with one row per fate per stage (including the dead fate)
and the following columns:

- from_stage:

  Character. The source stage (column of the matrix).

- to_stage:

  Character. The destination stage, including `"dead"`.

- mean:

  Numeric. Posterior mean transition probability.

- lower:

  Numeric. Lower bound of the credible interval.

- upper:

  Numeric. Upper bound of the credible interval.

## See also

[`plot_transition_CrI`](https://atiretoo.github.io/raretrans/reference/plot_transition_CrI.md)
for visualising the output.

## Examples

``` r
T_mat <- matrix(c(0.5, 0.3, 0.0,
                  0.2, 0.4, 0.1,
                  0.0, 0.1, 0.7), nrow = 3, ncol = 3)
F_mat <- matrix(c(0.0, 0.0, 1.5,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0), nrow = 3, ncol = 3)
TF <- list(T = T_mat, F = F_mat)
N  <- c(10, 5, 8)

# Default 95% credible intervals
cri <- transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"))
cri
#>    from_stage to_stage       mean        lower     upper
#> 1    plantula plantula 0.47727273 2.050069e-01 0.7573488
#> 2    plantula juvenile 0.29545455 7.933447e-02 0.5813692
#> 3    plantula    adult 0.02272727 2.540382e-08 0.1523734
#> 4    plantula     dead 0.20454545 3.411207e-02 0.4745221
#> 5    juvenile plantula 0.20833333 1.210489e-02 0.5781745
#> 6    juvenile juvenile 0.37500000 7.220923e-02 0.7550201
#> 7    juvenile    adult 0.12500000 1.277616e-03 0.4572469
#> 8    juvenile     dead 0.29166667 3.612337e-02 0.6744512
#> 9       adult plantula 0.02777778 3.146432e-08 0.1851079
#> 10      adult juvenile 0.11666667 3.866234e-03 0.3785689
#> 11      adult    adult 0.65000000 3.323634e-01 0.9050806
#> 12      adult     dead 0.20555556 2.575628e-02 0.5055951

# 90% credible intervals
transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"), ci = 0.90)
#>    from_stage to_stage       mean        lower     upper
#> 1    plantula plantula 0.47727273 2.419147e-01 0.7175380
#> 2    plantula juvenile 0.29545455 1.019478e-01 0.5328783
#> 3    plantula    adult 0.02272727 4.064623e-07 0.1100185
#> 4    plantula     dead 0.20454545 4.799268e-02 0.4238845
#> 5    juvenile plantula 0.20833333 2.140606e-02 0.5095062
#> 6    juvenile juvenile 0.37500000 1.007898e-01 0.6997954
#> 7    juvenile    adult 0.12500000 3.234720e-03 0.3843419
#> 8    juvenile     dead 0.29166667 5.490716e-02 0.6119090
#> 9       adult plantula 0.02777778 5.034306e-07 0.1344012
#> 10      adult juvenile 0.11666667 7.575849e-03 0.3214247
#> 11      adult    adult 0.65000000 3.826318e-01 0.8776073
#> 12      adult     dead 0.20555556 3.851405e-02 0.4493039
```
