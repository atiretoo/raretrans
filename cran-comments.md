## Release summary

This is a minor release. The main change is a bug fix in `fill_fertility()`
where vector recycling caused the β posterior to be miscomputed when
`priorweight > 0`. See NEWS.md for the full list of changes.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: macOS Tahoe 26.4.1 (aarch64-apple-darwin20), R 4.5.3
* GitHub Actions:
  - macos-latest (R release)
  - windows-latest (R release)
  - ubuntu-latest (R release, R-devel, R-oldrel-1)
* win-builder (planned: R-devel via `devtools::check_win_devel()`)

## Downstream dependencies

There are no reverse dependencies for this package on CRAN.
