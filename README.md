
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raretrans

[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/atiretoo/raretrans/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/atiretoo/raretrans/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/atiretoo/raretrans/branch/master/graph/badge.svg)](https://codecov.io/gh/atiretoo/raretrans)

Functions to create matrix population models from a combination of data
on stage/age transitions and Bayesian prior information. This
compensates for structural problems caused by missing observations of
rare transitions (“holey matrices”).

Based on methods described in:

> Tremblay, R. L., Tyre, A. J., Pérez, M.-E., & Ackerman, J. D. (2021).
> Population projections from holey matrices: Using prior information to
> estimate rare transition events. *Ecological Modelling*, **447**,
> 109526. <https://doi.org/10.1016/j.ecolmodel.2021.109526>

## Requirements

- R \>= 4.1.0
- Core package dependencies: `ggplot2`, `rlang`
- Vignette dependencies: `tidyverse`, `popbio`, `huxtable`, `popdemo`,
  `googledrive`

> **Note:** The core functions work on any R \>= 4.1.0. Vignettes
> require R \>= 4.1.0 due to dependencies in `popdemo`, `dplyr`, and
> `purrr`. If you are on an older version of R, install without
> vignettes (see below) and read them online at
> <https://atiretoo.github.io/raretrans>.

## Installation

`raretrans` is not currently available from CRAN. Install from GitHub
with:

``` r
# install.packages("remotes")  # if needed
remotes::install_github("atiretoo/raretrans")
```

To install with vignettes (requires R \>= 4.1.0):

``` r
remotes::install_github("atiretoo/raretrans", build_vignettes = TRUE)
```

If vignette installation fails, install without them and read online:

``` r
remotes::install_github("atiretoo/raretrans", build_vignettes = FALSE)
```

The code and data used to produce the published paper are tagged
`v1.0.0`:

``` r
remotes::install_github("atiretoo/raretrans@v1.0.0")
```

## Overview

The main functions in this package are:

| Function | Purpose |
|----|----|
| `fill_transitions()` | Fill a transition matrix using a Dirichlet-multinomial prior |
| `fill_fertility()` | Fill a fertility matrix using a Gamma prior |
| `sim_transitions()` | Simulate posterior transition matrices |
| `transition_CrI()` | Compute credible intervals for all transition probabilities |
| `plot_transition_CrI()` | Visualise posterior means and credible intervals |
| `plot_transition_density()` | Visualise full posterior beta densities as a matrix plot |

All functions take a `TF` list of two matrices as their first argument:
a transition matrix `T` and a fertility matrix `F`. These can be
constructed by hand or created with
`popbio::projection.matrix(..., TF = TRUE)` from individual-level
transition data.

## Quick Start

``` r
library(raretrans)

# Example with *Lepanthes eltoroensis* data included in the package
data(L_elto)

# Build a TF list from individual transition data
TF <- popbio::projection.matrix(L_elto, TF = TRUE)

# Fill structural zeros and perfect (1.0) survival, 
# transitions, stasis and fertilities using uninformative priors
fill_transitions(TF)
fill_fertility(TF)
```

## Vignettes

| \# | Vignette | Topic |
|----|----|----|
| 01 | Quick start | Minimal worked example with *Lepanthes eltoroensis* |
| 02 | Introduction | Full workflow walkthrough |
| 03 | Single population | Working with one population and transition period |
| 04 | Credible intervals | `transition_CrI()` and plotting with *Cypripedium calceolus* |
| 05 | Prior information | Effect of priors on transitions and fertility |
| 06 | Animal matrices | Examples with animal population data |

Browse all vignettes at <https://atiretoo.github.io/raretrans>.

## Getting Help

- Report bugs or request features at
  <https://github.com/atiretoo/raretrans/issues>
- Read the full function reference at
  <https://atiretoo.github.io/raretrans/reference>

## Reference

Caswell, H. (2001). *Matrix Population Models: Construction, Analysis,
and Interpretation* (2nd ed.). Sinauer Associates.

## Code of Conduct

Please note that the `raretrans` project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.
