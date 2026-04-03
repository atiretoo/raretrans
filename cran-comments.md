## Test environments

* local: macOS Tahoe 26.4, R 4.5.2 (aarch64-apple-darwin20)
* GitHub Actions: ubuntu-latest (release, devel), windows-latest (release), macOS-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission notes

This is a resubmission of version 1.0.2. Changes since 1.0.1:

* Added three new exported functions: `transition_CrI()`, `plot_transition_CrI()`,
  and `plot_transition_density()`
* Replaced `tidyverse` in Suggests with the specific packages used
  (`dplyr`, `tidyr`, `purrr`, `tibble`)
* Removed `devtools` and `googledrive` from Suggests
* Cleaned up vignette build artefacts and developer-only code chunks
* All vignettes build cleanly with 0 errors, 0 warnings, 0 notes

## Downstream dependencies

There are currently no downstream dependencies for this package.
