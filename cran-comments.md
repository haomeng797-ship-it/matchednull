## Test environments

* local: macOS (Apple Silicon), R 4.5.2
* win-builder: R-devel, R-release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

* This is a minor release. `copula_null()` and `matched_null_test()` gain a
  `copula = "t"` option with a `df` argument, so that an apparent excess of
  clusters can be checked against a heavier-tailed null before it is read as
  evidence of cluster structure.
* The mclust, knitr, and rmarkdown packages are used only in Suggests
  (examples, tests, and the vignette are conditional on their availability).
* Reverse dependencies: none.
