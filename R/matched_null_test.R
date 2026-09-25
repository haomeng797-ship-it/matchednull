#' Test a cluster count against matched nulls
#'
#' Runs the user's own clustering pipeline on the real data and on `R`
#' matched-null twins of it (see [copula_null()]), and asks whether the real
#' result exceeds what the twins produce. The twins are a reference model
#' fitted to the data: every twin column is a permutation of the observed
#' values, so the margins are preserved exactly; the correlations are
#' reproduced approximately (see [copula_null()] for why); and the twins
#' contain no cluster structure by construction. The pipeline is supplied as
#' a function, so any procedure that returns a number, a BIC-selected
#' mixture, a k-means heuristic, a published typology's own workflow, can be
#' tested unchanged.
#'
#' `cluster_fn` must return one finite number on every call. If it returns
#' `NA`, `NaN`, `Inf`, a non-numeric value, or a vector of the wrong length,
#' on the real data or on any twin, the test stops with an error naming the
#' offending call. Invalid twins are never dropped and `R` is never quietly
#' reduced, because an exceedance computed from a subset of the twins would
#' not be the one that was asked for. Errors raised inside `cluster_fn`
#' propagate unchanged.
#'
#' `p_exceed` is a plug-in Monte Carlo exceedance quantity: one plus the
#' number of twins whose statistic reaches the real value, divided by
#' `R + 1`, under a reference fitted to the same data. It is not an exact
#' finite-sample p-value. The reference's correlations are estimated from
#' `x` and reproduced only approximately, so `p_exceed` is a comparison
#' against that fitted reference, not a test with guaranteed size.
#'
#' Reading the outcome. A real statistic above the twins' interval means the
#' data depart from this reference model, one that carries their margins,
#' approximately their correlations, and no cluster structure. It does not by
#' itself establish that clusters exist: structure the reference does not
#' carry, heavy-tailed dependence included, can produce the same departure,
#' and whether a particular departure registers depends on the statistic, the
#' pipeline and the sample size. A real statistic inside the interval means no
#' detected excess under this reference; it is not evidence that clusters are
#' absent. A statistic below the interval is reported as a diagnostic (the
#' pipeline finds less on the data than on its twins), not as a verdict.
#'
#' `copula = "t"` is a sensitivity variant. It gives the twins heavier-tailed
#' dependence at a `df` fixed by the caller (8, then 3, is a reasonable
#' pair), so a rerun shows whether a departure from the Gaussian twins
#' persists once tail co-movement is built into the reference. A departure
#' that persists is harder to attribute to tails; one that disappears is
#' consistent with tail dependence. Neither outcome by itself establishes or
#' rules out clusters. The Gaussian run is the primary analysis and the t run
#' is reported alongside it, not in its place.
#'
#' @param x A numeric matrix or data frame, complete cases only.
#' @param cluster_fn A function taking a data matrix and returning one finite
#'   number: the statistic to be tested, typically the selected number of
#'   clusters (but any scalar summary of clustering strength works). A return
#'   value that is missing, non-finite, non-numeric, or not of length one is
#'   an error, on the real data and on every twin.
#' @param R Number of matched-null twins to draw (default 200). The result
#'   always holds exactly `R` twin statistics.
#' @param probs Lower and upper quantiles of the twin statistics that define
#'   the reference interval (default `c(.025, .975)`).
#' @param copula,df Dependence family of the twins and t degrees of freedom,
#'   passed to [copula_null()]. Defaults to `"gaussian"`, which is the
#'   reference model of the method; `copula = "t"` is a sensitivity variant
#'   with heavier-tailed dependence, run in addition to the Gaussian test
#'   rather than in place of it, and `df` is fixed by the caller, never
#'   estimated from the data.
#' @param ridge Passed to [copula_null()].
#' @param parallel Evaluate the twin loop through \pkg{future.apply}, so that the
#'   backend is whatever `future::plan()` the caller has set. Because the twins
#'   are independent and the null construction itself is a negligible share of
#'   the cost, the speed-up is close to the cost of the user's own pipeline
#'   divided by the number of workers. Seeded results are reproducible and do
#'   not depend on how many workers run them. They do, however, differ from the
#'   `parallel = FALSE` stream, which is left byte-identical to earlier versions,
#'   so a seeded analysis should not switch settings partway through.
#'
#' @return An object of class `"matched_null_test"`: a list with the real
#'   statistic (`real`), the `R` twin statistics (`null`), the reference
#'   interval (`interval`, the `probs` quantiles of `null`), the plug-in
#'   Monte Carlo exceedance quantity (`p_exceed`, computed as
#'   `(1 + sum(null >= real)) / (R + 1)`; not an exact finite-sample p-value,
#'   see Details), whether the real statistic lies inside the interval
#'   (`within`), the requested `R`, and the reference family used (`copula`,
#'   with `df` when it is `"t"`). Reproducibility is the caller's: set a seed
#'   before calling.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mclust", quietly = TRUE)) {
#'   suppressPackageStartupMessages(library(mclust))
#'   # a pipeline: how many components does BIC select?
#'   pick_k <- function(d) Mclust(d, G = 1:5, verbose = FALSE)$G
#'   set.seed(42)
#'   x <- matrix(rnorm(500 * 4), 500, 4)  # typeless data
#'   matched_null_test(x, pick_k, R = 30)
#'
#'   # sensitivity variant with heavier-tailed dependence in the reference:
#'   # matched_null_test(x, pick_k, R = 30, copula = "t", df = 8)
#' }
#' }
#'
#' @export
matched_null_test <- function(x, cluster_fn, R = 200, probs = c(.025, .975),
                              copula = c("gaussian", "t"), df = 8,
                              ridge = 1e-6, parallel = FALSE) {
  x <- as.matrix(x)
  copula <- match.arg(copula)
  if (!is.function(cluster_fn)) stop("`cluster_fn` must be a function.", call. = FALSE)
  if (!is.numeric(R) || length(R) != 1L || R < 1) stop("`R` must be a positive number.", call. = FALSE)
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel))
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)

  check <- .check_statistic   # bound here so the twin closure carries it to parallel workers
  real <- check(cluster_fn(x), "the real data")

  draw_one <- function(i) {
    check(
      cluster_fn(copula_null(x, copula = copula, df = df, ridge = ridge)),
      sprintf("null twin %d of %d", as.integer(i), as.integer(R))
    )
  }

  nulls <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("`parallel = TRUE` requires the future.apply package.", call. = FALSE)
    future.apply::future_vapply(seq_len(R), draw_one, numeric(1), future.seed = TRUE)
  } else {
    vapply(seq_len(R), draw_one, numeric(1))
  }

  q <- stats::quantile(nulls, probs)
  out <- list(
    real     = real,
    null     = nulls,
    interval = q,
    p_exceed = (1 + sum(nulls >= real)) / (R + 1),
    within   = real >= q[1] && real <= q[2],
    R        = R,
    copula   = copula,
    df       = if (copula == "t") df else NULL
  )
  class(out) <- "matched_null_test"
  out
}

# Validate the raw return value of `cluster_fn` before any coercion: it must
# be one finite number. `what` names the call ("the real data" or
# "null twin r of R") so that the error points at the run that failed.
# Nothing invalid is ever dropped and carried on from; the test stops.
.check_statistic <- function(value, what) {
  problem <- if (!is.numeric(value)) {
    sprintf("returned an object of class <%s> and length %d",
            paste(class(value), collapse = "/"), length(value))
  } else if (length(value) != 1L) {
    sprintf("returned a numeric vector of length %d", length(value))
  } else if (!is.finite(value)) {
    sprintf("returned %s", format(value))
  } else {
    NULL
  }
  if (!is.null(problem))
    stop(sprintf("`cluster_fn` %s on %s; it must return one finite number.",
                 problem, what), call. = FALSE)
  unname(as.numeric(value))
}

#' @export
print.matched_null_test <- function(x, ...) {
  fam <- if (is.null(x$copula) || x$copula == "gaussian") "Gaussian"
         else paste0("t (df = ", x$df, ")")
  cat("Matched-null test (", x$R, " ", fam, " null twins)\n", sep = "")
  cat("  real statistic:               ", x$real, "\n", sep = "")
  cat("  null interval:                [", paste(round(x$interval, 2), collapse = ", "),
      "]\n", sep = "")
  # the numerator counts twins that reach the real value, so the label names
  # the direction actually counted: null >= real, not real >= null
  cat("  MC exceedance (null >= real): ", round(x$p_exceed, 3), "\n", sep = "")
  cat("  outcome:                      ",
      if (x$within) "no detected excess (within the reference interval)"
      else if (x$real > x$interval[2]) "departure from the reference model (above the interval; kind not settled by this test)"
      else "below the reference range (a diagnostic, not a verdict)", "\n", sep = "")
  invisible(x)
}
