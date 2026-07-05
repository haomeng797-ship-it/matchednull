#' Draw one Gaussian-copula matched null ("null twin") of a dataset
#'
#' Builds a synthetic twin of `x` that preserves every marginal distribution
#' exactly and the correlation matrix to within sampling error, while
#' containing no cluster structure by construction: all dependence in the twin
#' is Gaussian. Clustering found in real data but not in its twins must come
#' from structure beyond the margins and covariance; clustering found equally
#' in both was never more than the data's shape.
#'
#' The construction is a rank reordering in the tradition of Iman and Conover
#' (1982): draw a Gaussian sample with the data's correlation matrix, then
#' replace each column with the sorted real values laid down in the rank order
#' of the Gaussian column. Every real value is reused exactly once per column,
#' which is why the margins match exactly; the rank correlation is matched
#' exactly and the Pearson correlation to within sampling error.
#'
#' @param x A numeric matrix or data frame (rows = observations, columns =
#'   variables), complete cases only, at least two rows and two columns.
#' @param ridge Small value added to the diagonal of the correlation matrix
#'   only if it is not positive definite (default `1e-6`).
#'
#' @return A numeric matrix of the same dimensions as `x`: one matched-null
#'   draw. Each column contains exactly the values of the corresponding column
#'   of `x`, rearranged.
#'
#' @references
#' Iman, R. L., & Conover, W. J. (1982). A distribution-free approach to
#' inducing rank correlation among input variables. *Communications in
#' Statistics - Simulation and Computation, 11*(3), 311-334.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200 * 3), 200, 3) %*% chol(matrix(c(1, .5, .3,
#'                                                        .5, 1, .4,
#'                                                        .3, .4, 1), 3, 3))
#' twin <- copula_null(x)
#' # margins identical:
#' all(sort(twin[, 1]) == sort(x[, 1]))
#' # correlations close:
#' round(cor(x) - cor(twin), 2)
#'
#' @export
copula_null <- function(x, ridge = 1e-6) {
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (anyNA(x)) stop("`x` contains missing values; supply complete cases.", call. = FALSE)
  n <- nrow(x); p <- ncol(x)
  if (n < 2L || p < 2L) stop("`x` needs at least 2 rows and 2 columns.", call. = FALSE)

  C <- stats::cor(x)
  L <- tryCatch(chol(C), error = function(e) chol(C + diag(ridge, p)))
  Z <- matrix(stats::rnorm(n * p), n, p) %*% L

  out <- matrix(0, n, p, dimnames = dimnames(x))
  for (j in seq_len(p)) {
    s <- sort(x[, j])
    out[, j] <- s[rank(Z[, j], ties.method = "first")]
  }
  out
}
