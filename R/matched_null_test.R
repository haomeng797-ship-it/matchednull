#' Test a cluster count against Gaussian-copula matched nulls
#'
#' Runs the user's own clustering pipeline on the real data and on `R`
#' matched-null twins of it (see [copula_null()]), and asks whether the real
#' result exceeds what the twins, which share the data's margins and
#' covariance but contain no cluster structure, produce. The pipeline is
#' supplied as a function, so any procedure that returns a number, a
#' BIC-selected mixture, a k-means heuristic, a published typology's own
#' workflow, can be tested unchanged.
#'
#' @param x A numeric matrix or data frame, complete cases only.
#' @param cluster_fn A function taking a data matrix and returning a single
#'   number: the statistic to be tested, typically the selected number of
#'   clusters (but any scalar summary of clustering strength works).
#' @param R Number of matched-null twins to draw (default 200).
#' @param probs Lower and upper quantiles of the null distribution used for
#'   the interval verdict (default `c(.025, .975)`).
#' @param ridge Passed to [copula_null()].
#'
#' @return An object of class `"matched_null_test"`: a list with the real
#'   statistic (`real`), the null draws (`null`), the null interval
#'   (`interval`), a one-sided Monte Carlo p-value for exceedance
#'   (`p_exceed`), and the interval verdict (`within`). Reproducibility is the
#'   caller's: set a seed before calling.
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
#' }
#' }
#'
#' @export
matched_null_test <- function(x, cluster_fn, R = 200, probs = c(.025, .975),
                              ridge = 1e-6) {
  x <- as.matrix(x)
  if (!is.function(cluster_fn)) stop("`cluster_fn` must be a function.", call. = FALSE)
  if (!is.numeric(R) || length(R) != 1L || R < 1) stop("`R` must be a positive number.", call. = FALSE)

  real <- as.numeric(cluster_fn(x))
  if (length(real) != 1L || is.na(real))
    stop("`cluster_fn` must return a single non-missing number.", call. = FALSE)

  nulls <- vapply(seq_len(R), function(i) {
    as.numeric(cluster_fn(copula_null(x, ridge = ridge)))
  }, numeric(1))

  q <- stats::quantile(nulls, probs, na.rm = TRUE)
  out <- list(
    real     = real,
    null     = nulls,
    interval = q,
    p_exceed = (1 + sum(nulls >= real, na.rm = TRUE)) / (R + 1),
    within   = real >= q[1] && real <= q[2],
    R        = R
  )
  class(out) <- "matched_null_test"
  out
}

#' @export
print.matched_null_test <- function(x, ...) {
  cat("Matched-null test (", x$R, " null twins)\n", sep = "")
  cat("  real statistic:      ", x$real, "\n", sep = "")
  cat("  null interval:       [", paste(round(x$interval, 2), collapse = ", "),
      "]\n", sep = "")
  cat("  p (real >= nulls):   ", round(x$p_exceed, 3), "\n", sep = "")
  cat("  verdict:             ",
      if (x$within) "null-like (within the twins' interval)"
      else if (x$real > x$interval[2]) "exceeds the null (beyond margins + covariance)"
      else "below the null interval", "\n", sep = "")
  invisible(x)
}
