## Test environments

* local: macOS (Apple Silicon), R 4.5.2
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  release, and oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

* This is a patch release. `matched_null_test()` now stops with an
  informative error when the user's clustering function returns an invalid
  value, instead of silently leaving that twin out of the exceedance count.
  The printed summary relabels the exceedance line, whose old label named the
  wrong direction, and renames the three outcomes. The documentation now says
  that the correlations are reproduced approximately rather than exactly, and
  the citation file gives the current title of the accompanying manuscript.
  Results for valid inputs are unchanged.
* The mclust, knitr, rmarkdown, future, and future.apply packages are used
  only in Suggests and are checked before use.
* Reverse dependencies: none.
