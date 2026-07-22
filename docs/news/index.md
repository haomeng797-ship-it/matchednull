# Changelog

## matchednull (development version)

- [`copula_null()`](https://haomeng797-ship-it.github.io/matchednull/reference/copula_null.md)
  and
  [`matched_null_test()`](https://haomeng797-ship-it.github.io/matchednull/reference/matched_null_test.md)
  gain `copula = "t"` and `df`: t-copula twins preserve the same margins
  and correlations but add tail dependence. A verdict of “exceeds the
  Gaussian null” can now be stress-tested against heavier-tailed
  dependence (`df = 8`, then `3`) before it is read as evidence of
  types.
- [`matched_null_test()`](https://haomeng797-ship-it.github.io/matchednull/reference/matched_null_test.md)
  results record the null family used (`$copula`, `$df`), and the
  printed summary names it.

## matchednull 0.1.0

CRAN release: 2026-07-21

- First release:
  [`copula_null()`](https://haomeng797-ship-it.github.io/matchednull/reference/copula_null.md)
  builds Gaussian-copula matched-null twins;
  [`matched_null_test()`](https://haomeng797-ship-it.github.io/matchednull/reference/matched_null_test.md)
  tests any clustering pipeline against them.
