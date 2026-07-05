test_that("margins are preserved exactly", {
  set.seed(1)
  x <- matrix(rnorm(100 * 3), 100, 3)
  twin <- copula_null(x)
  for (j in 1:3) expect_identical(sort(twin[, j]), sort(x[, j]))
})

test_that("dimensions and dimnames are preserved", {
  set.seed(2)
  x <- matrix(rnorm(60 * 2), 60, 2, dimnames = list(NULL, c("a", "b")))
  twin <- copula_null(x)
  expect_identical(dim(twin), dim(x))
  expect_identical(colnames(twin), c("a", "b"))
})

test_that("correlations are close to the original", {
  set.seed(3)
  x <- matrix(rnorm(2000 * 3), 2000, 3) %*%
    chol(matrix(c(1, .6, .3, .6, 1, .4, .3, .4, 1), 3, 3))
  twin <- copula_null(x)
  expect_lt(max(abs(cor(x) - cor(twin))), 0.06)
})

test_that("bad input fails loudly", {
  expect_error(copula_null(matrix(c(1, NA, 3, 4), 2, 2)), "missing")
  expect_error(copula_null(matrix(letters[1:4], 2, 2)), "numeric")
  expect_error(copula_null(matrix(1:5, 5, 1)), "2 rows and 2 columns")
})

test_that("a data.frame input works", {
  set.seed(4)
  d <- data.frame(u = rnorm(50), v = rnorm(50))
  expect_silent(twin <- copula_null(d))
  expect_identical(sort(twin[, "u"]), sort(d$u))
})
