# Changelog

## raretrans 1.0.2

### Maintainer

- Raymond Tremblay (University of Puerto Rico at Humacao) is now the
  package maintainer. Andrew Tyre remains as author.

### New functions

- [`transition_CrI()`](https://atiretoo.github.io/raretrans/reference/transition_CrI.md):
  computes marginal posterior beta credible intervals for all entries of
  the transition matrix, returning a tidy data frame.
- [`plot_transition_CrI()`](https://atiretoo.github.io/raretrans/reference/plot_transition_CrI.md):
  visualises posterior mean transition probabilities and credible
  intervals as a point-range plot, with optional exclusion of the dead
  fate (`include_dead`).
- [`plot_transition_density()`](https://atiretoo.github.io/raretrans/reference/plot_transition_density.md):
  visualises the full posterior beta density for each matrix entry,
  arranged as an n x n matrix plot mirroring the structure of the
  population projection matrix, with optional exclusion of the dead fate
  (`include_dead`).

### Bug fixes

- Fixed
  [`fill_transitions()`](https://atiretoo.github.io/raretrans/reference/fill_transitions.md)
  and
  [`sim_transitions()`](https://atiretoo.github.io/raretrans/reference/sim_transitions.md):
  replaced `missing(P)` with `is.null(P)` so that passing `P = NULL`
  correctly triggers the uniform prior default.

### Vignettes

- Added `introduction` vignette: a beginner-friendly walkthrough of the
  full `raretrans` workflow using a published *Chamaedorea elegans*
  matrix (Burns et al. 2013 / COMPADRE). Covers entering matrices
  directly, uniform and informative priors,
  [`fill_transitions()`](https://atiretoo.github.io/raretrans/reference/fill_transitions.md),
  [`fill_fertility()`](https://atiretoo.github.io/raretrans/reference/fill_fertility.md),
  [`sim_transitions()`](https://atiretoo.github.io/raretrans/reference/sim_transitions.md),
  [`transition_CrI()`](https://atiretoo.github.io/raretrans/reference/transition_CrI.md),
  [`plot_transition_CrI()`](https://atiretoo.github.io/raretrans/reference/plot_transition_CrI.md),
  and
  [`plot_transition_density()`](https://atiretoo.github.io/raretrans/reference/plot_transition_density.md).
- Fixed `onepopperiod` vignette: removed
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
  call that would fail during CRAN’s vignette build.

### Documentation

- Completed full documentation for the `L_elto` dataset, including
  accurate descriptions of all 13 variables, stage code definitions, and
  references.
- Added `@examples` to all exported functions.
- Fixed typo in package-level documentation (“observatons” -\>
  “observations”).

### DESCRIPTION

- Fixed duplicate `URL` and `BugReports` fields.
- Updated URLs to point to the current repository (`atiretoo`).
- Replaced `tidyverse` in `Suggests` with the specific packages used
  (`dplyr`, `purrr`, `tibble`).
- Added `ggplot2` and `rlang` to `Imports`.
- Version bumped from 1.0.1.9000 to 1.0.2.

## raretrans 1.0.1.9000

### Breaking changes

- Fixed Issue [\#2](https://github.com/atiretoo/raretrans/issues/2),
  which means code that worked with v1.0.0 will not work without
  warnings, and will likely produce incorrect results.

### Additions

- Added tests for all major functions

### Minor fixes

- Added a `NEWS.md` file to track changes to the package.
- Used tidy styling, added issue template, contribution guidelines etc.

## raretrans 1.0.0

- First version of package, used to create figures in Tremblay et
  al. (submitted)
