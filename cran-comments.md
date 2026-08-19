## Test environments

* local: macOS (Apple Silicon), R 4.5.2
* win-builder: R-devel, R-release

## R CMD check results

0 errors | 0 warnings | 1 note (the expected new-maintainer note, explained below)

## Comments

* The maintainer email changes with this release, from meng10@upenn.edu to
  haomeng797@gmail.com. Both addresses are mine: I am moving the package to my
  permanent personal address while the university one, which I will lose after
  graduating, still works and can be used to confirm the change.
* This is a minor release. `copula_null()` and `matched_null_test()` gain a
  `copula = "t"` option with a `df` argument, so that an apparent excess of
  clusters can be checked against a heavier-tailed null before it is read as
  evidence of cluster structure. `matched_null_test()` also gains optional
  parallel execution through 'future.apply'.
* The mclust, knitr, rmarkdown, future, and future.apply packages are used
  only in Suggests and are checked before use.
* Reverse dependencies: none.
