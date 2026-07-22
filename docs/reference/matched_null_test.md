# Test a cluster count against matched nulls

Runs the user's own clustering pipeline on the real data and on `R`
matched-null twins of it (see
[`copula_null()`](https://haomeng797-ship-it.github.io/matchednull/reference/copula_null.md)),
and asks whether the real result exceeds what the twins, which share the
data's margins and covariance but contain no cluster structure, produce.
The pipeline is supplied as a function, so any procedure that returns a
number, a BIC-selected mixture, a k-means heuristic, a published
typology's own workflow, can be tested unchanged.

## Usage

``` r
matched_null_test(
  x,
  cluster_fn,
  R = 200,
  probs = c(0.025, 0.975),
  copula = c("gaussian", "t"),
  df = 8,
  ridge = 1e-06
)
```

## Arguments

- x:

  A numeric matrix or data frame, complete cases only.

- cluster_fn:

  A function taking a data matrix and returning a single number: the
  statistic to be tested, typically the selected number of clusters (but
  any scalar summary of clustering strength works).

- R:

  Number of matched-null twins to draw (default 200).

- probs:

  Lower and upper quantiles of the null distribution used for the
  interval verdict (default `c(.025, .975)`).

- copula, df:

  Dependence family of the twins and t degrees of freedom, passed to
  [`copula_null()`](https://haomeng797-ship-it.github.io/matchednull/reference/copula_null.md).
  Defaults to `"gaussian"`.

- ridge:

  Passed to
  [`copula_null()`](https://haomeng797-ship-it.github.io/matchednull/reference/copula_null.md).

## Value

An object of class `"matched_null_test"`: a list with the real statistic
(`real`), the null draws (`null`), the null interval (`interval`), a
one-sided Monte Carlo p-value for exceedance (`p_exceed`), the interval
verdict (`within`), and the null family used (`copula`, with `df` when
it is `"t"`). Reproducibility is the caller's: set a seed before
calling.

## Details

A verdict of "exceeds the null" under the default Gaussian twins
licenses only "structure beyond margins and correlations", it does not
by itself license types: heavy-tailed dependence also exceeds a Gaussian
null. To separate the two, rerun with `copula = "t"` (df 8, then 3). A
result that survives the t twins as well is harder to attribute to
tails; a result the t twins reproduce was tail dependence, not types.

## Examples

``` r
# \donttest{
if (requireNamespace("mclust", quietly = TRUE)) {
  suppressPackageStartupMessages(library(mclust))
  # a pipeline: how many components does BIC select?
  pick_k <- function(d) Mclust(d, G = 1:5, verbose = FALSE)$G
  set.seed(42)
  x <- matrix(rnorm(500 * 4), 500, 4)  # typeless data
  matched_null_test(x, pick_k, R = 30)

  # stress test of an exceedance verdict against heavier tails:
  # matched_null_test(x, pick_k, R = 30, copula = "t", df = 8)
}
#> Matched-null test (30 Gaussian null twins)
#>   real statistic:      1
#>   null interval:       [1, 1]
#>   p (real >= nulls):   1
#>   verdict:             null-like (within the twins' interval)
# }
```
