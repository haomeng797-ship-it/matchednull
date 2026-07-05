test_that("mechanics: structure, p-value bounds, printing", {
  set.seed(10)
  x <- matrix(rnorm(80 * 2), 80, 2)
  fixed <- function(d) 3            # a pipeline that always answers 3
  r <- matched_null_test(x, fixed, R = 19)
  expect_s3_class(r, "matched_null_test")
  expect_identical(r$real, 3)
  expect_length(r$null, 19)
  expect_true(r$p_exceed > 0 && r$p_exceed <= 1)
  expect_true(r$within)             # real == every null, so inside the interval
  expect_output(print(r), "null-like")
})

test_that("bad input fails loudly", {
  x <- matrix(rnorm(40), 20, 2)
  expect_error(matched_null_test(x, cluster_fn = 5), "must be a function")
  expect_error(matched_null_test(x, function(d) c(1, 2), R = 5), "single non-missing")
  expect_error(matched_null_test(x, function(d) 1, R = 0), "positive")
})

test_that("behaviour: quiet on typeless data, fires on dependence types", {
  skip_if_not_installed("mclust")
  suppressPackageStartupMessages(require(mclust, quietly = TRUE))
  pick_k <- function(d)
    mclust::Mclust(d, G = 1:4, modelNames = "VVV", verbose = FALSE)$G

  set.seed(42)
  x0 <- matrix(rnorm(400 * 4), 400, 4) %*% chol(diag(4) * .5 + .5)
  set.seed(7)
  quiet <- matched_null_test(x0, pick_k, R = 19)
  expect_true(quiet$within)

  set.seed(42)
  z <- sample(2, 400, replace = TRUE); X <- matrix(rnorm(400 * 4), 400, 4)
  L1 <- chol(matrix(c(1, .85, .85, 1), 2, 2))
  L2 <- chol(matrix(c(1, -.85, -.85, 1), 2, 2))
  X[z == 1, 1:2] <- X[z == 1, 1:2] %*% L1; X[z == 1, 3:4] <- X[z == 1, 3:4] %*% L1
  X[z == 2, 1:2] <- X[z == 2, 1:2] %*% L2; X[z == 2, 3:4] <- X[z == 2, 3:4] %*% L2
  set.seed(7)
  fired <- matched_null_test(X, pick_k, R = 19)
  expect_gt(fired$real, fired$interval[2])
})
