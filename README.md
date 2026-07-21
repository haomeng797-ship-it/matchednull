# matchednull

Clustering methods will always give you clusters, even when the data are one
smooth cloud. That makes results hard to trust: is the k your pipeline reports
a property of the data, or just of the method?

matchednull helps you check. It builds a twin of your dataset with the same
margins and the same correlations but no cluster structure, then runs your own
pipeline on both. If the real data give you more clusters than the twins do,
that's real signal. If not, now you know before a reviewer asks.

`copula_null()` builds the twin: every marginal distribution is preserved
exactly and the correlation matrix to within sampling error, with no cluster
structure by construction. `matched_null_test()` runs any clustering pipeline,
supplied as a function, on the real data and on `R` twins, and asks whether
the real result stands out.

## Installation

```r
install.packages("matchednull")

# dev version, has the t-copula stress test
remotes::install_github("haomeng797-ship-it/matchednull")
```

## Quick example

```r
library(matchednull)
library(mclust)

pick_k <- function(d) Mclust(d, G = 1:5, verbose = FALSE)$G

set.seed(1)
x <- matrix(rnorm(500 * 4), 500, 4)   # typeless data
set.seed(2)
matched_null_test(x, pick_k, R = 50)
#> Matched-null test (50 Gaussian null twins)
#>   real statistic:      1
#>   null interval:       [1, 1]
#>   p (real >= nulls):   1
#>   verdict:             null-like (within the twins' interval)
```

See `vignette("matchednull")` for the full walk-through, including a case
where genuine types hide in the dependence structure and the test fires.

In our calibration runs on skewed but clusterless data, standard criteria
(BIC, the bootstrap LRT) reported structure in every dataset. The matched-null
test flagged none of them.

## Reference

The method, its positive controls, and its false-positive calibration are
described in the accompanying paper *Types Without Taxa* (Meng, 2026;
preregistration: <https://osf.io/2ekcg>).
