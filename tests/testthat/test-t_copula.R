test_that("t twins preserve margins exactly and correlations closely", {
  set.seed(21)
  x <- matrix(rnorm(2000 * 3), 2000, 3) %*%
    chol(matrix(c(1, .6, .3, .6, 1, .4, .3, .4, 1), 3, 3))
  twin <- copula_null(x, copula = "t", df = 5)
  for (j in 1:3) expect_identical(sort(twin[, j]), sort(x[, j]))
  expect_lt(max(abs(cor(x) - cor(twin))), 0.08)
})

test_that("bad df fails loudly, and df is ignored for gaussian", {
  x <- matrix(rnorm(40), 20, 2)
  expect_error(copula_null(x, copula = "t", df = -1), "positive")
  expect_error(copula_null(x, copula = "t", df = c(3, 8)), "single")
  expect_error(copula_null(x, copula = "t", df = NA), "positive|single")
  expect_silent(copula_null(x, df = -99))  # gaussian path never touches df
})

test_that("t twins carry tail dependence that gaussian twins lack", {
  set.seed(22)
  n <- 4000
  Z <- matrix(rnorm(n * 2), n, 2) %*% chol(matrix(c(1, .5, .5, 1), 2, 2))
  x <- Z / sqrt(rchisq(n, 3) / 3)  # t-copula data: one population, no clusters
  joint_tail <- function(d) {
    mean(d[, 1] > quantile(d[, 1], .95) & d[, 2] > quantile(d[, 2], .95))
  }
  real <- joint_tail(x)
  set.seed(23)
  g  <- replicate(40, joint_tail(copula_null(x)))
  set.seed(24)
  t3 <- replicate(40, joint_tail(copula_null(x, copula = "t", df = 3)))
  # the stress-test logic in miniature: gaussian twins cannot reproduce the
  # joint-tail mass, t twins can
  expect_gt(real, max(g))
  expect_gt(mean(t3), mean(g) * 1.5)
  expect_lt(abs(mean(t3) - real), real - mean(g))
})

test_that("matched_null_test passes the copula through and records it", {
  set.seed(25)
  x <- matrix(rnorm(300 * 2), 300, 2)
  joint_tail <- function(d) {
    mean(d[, 1] > quantile(d[, 1], .9) & d[, 2] > quantile(d[, 2], .9))
  }
  r <- matched_null_test(x, joint_tail, R = 19, copula = "t", df = 4)
  expect_identical(r$copula, "t")
  expect_identical(r$df, 4)
  expect_output(print(r), "t \\(df = 4\\)")
  g <- matched_null_test(x, joint_tail, R = 19)
  expect_identical(g$copula, "gaussian")
  expect_null(g$df)
  expect_output(print(g), "Gaussian")
})
