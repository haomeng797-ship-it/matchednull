# matchednull

Matched-null tests for cluster-count claims: does a reported number of
clusters or "types" exceed what the data's own margins and covariance already
produce?

`copula_null()` builds a synthetic twin of a dataset that preserves every
marginal distribution exactly and the correlation matrix to within sampling
error, while containing no cluster structure by construction.
`matched_null_test()` runs any clustering pipeline, supplied as a function, on
the real data and on `R` twins, and asks whether the real result stands out.

## Installation

```r
# development version
remotes::install_github("haomeng797-ship-it/matchednull")
```

A CRAN release is planned.

## Quick example

```r
library(matchednull)
library(mclust)

pick_k <- function(d) Mclust(d, G = 1:5, verbose = FALSE)$G

set.seed(1)
x <- matrix(rnorm(500 * 4), 500, 4)   # typeless data
set.seed(2)
matched_null_test(x, pick_k, R = 50)  # verdict: null-like
```

See `vignette("matchednull")` for the full walk-through, including a case
where genuine types hide in the dependence structure and the test fires.

## Reference

The method, its positive controls, and its false-positive calibration are
described in the accompanying paper *Types Without Taxa* (Meng, 2026;
preregistration: <https://osf.io/2ekcg>).
