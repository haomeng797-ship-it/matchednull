# matchednull 0.2.2

* `matched_null_test()` now stops with an informative error when `cluster_fn`
  returns an invalid value (missing, non-finite, non-numeric, or not of length
  one), on the real data or on any twin. The error names the call ("the real
  data" or "null twin r of R"). Previously an invalid twin was silently dropped
  from the numerator of `p_exceed` while the denominator kept the requested
  `R`, and the printed summary still reported all `R` twins, so a run in which
  most twins had failed could print a p-value and an outcome as if nothing had
  gone wrong. Errors raised inside `cluster_fn` propagate unchanged, as before.
  Numerical results for valid inputs are unchanged: seeded runs reproduce
  exactly.
* The printed exceedance line is relabelled. It read "p (real >= nulls)", but
  the numerator counts twins whose statistic reaches the real value, so the
  label named the wrong direction. It now reads
  "MC exceedance (null >= real)", and the other printed labels are re-padded
  to keep the block aligned. The value is unchanged:
  `(1 + sum(null >= real)) / (R + 1)`, as before.
* The printed outcome line is renamed. The three outcomes now read "no
  detected excess (within the reference interval)", "departure from the
  reference model (above the interval; kind not settled by this test)", and
  "below the reference range (a diagnostic, not a verdict)". They replace
  "null-like (within the twins' interval)", "exceeds the null (beyond margins
  + covariance)", and "below the null interval". Returned values are
  unchanged.
* Documentation corrected. The twins preserve the margins exactly (each column
  is a permutation of the observed values) but reproduce the correlations only
  approximately: the observed correlation matrix is used as the latent
  correlation of the copula, and non-normal margins change the output Pearson
  correlations, so the earlier "to within sampling error" wording was wrong.
  `p_exceed` is now described as a plug-in Monte Carlo exceedance quantity
  under a reference fitted to the same data, not an exact finite-sample
  p-value, and `copula = "t"` as a sensitivity variant with heavier-tailed
  dependence that does not by itself establish or rule out clusters.

# matchednull 0.2.1

* `matched_null_test()` gains `parallel`. With `parallel = TRUE` the twin
  loop runs through 'future.apply', so the backend is whatever `future::plan()`
  the caller has set. Seeded results are reproducible and independent of the
  number of workers. The default `parallel = FALSE` path is unchanged, so
  existing seeded analyses reproduce exactly.
* The maintainer address now uses a permanent personal inbox.

# matchednull 0.2.0

* `copula_null()` and `matched_null_test()` gain `copula = "t"` and `df`:
  t-copula twins preserve the same margins and correlations but add tail
  dependence. A verdict of "exceeds the Gaussian null" can now be
  stress-tested against heavier-tailed dependence (`df = 8`, then `3`)
  before it is read as evidence of types.
* `matched_null_test()` results record the null family used (`$copula`,
  `$df`), and the printed summary names it.

# matchednull 0.1.0

* First release: `copula_null()` builds Gaussian-copula matched-null twins;
  `matched_null_test()` tests any clustering pipeline against them.
