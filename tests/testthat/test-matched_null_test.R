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
  expect_output(print(r), "no detected excess")
})

test_that("bad input fails loudly", {
  x <- matrix(rnorm(40), 20, 2)
  expect_error(matched_null_test(x, cluster_fn = 5), "must be a function")
  expect_error(matched_null_test(x, function(d) c(1, 2), R = 5), "finite number")
  expect_error(matched_null_test(x, function(d) 1, R = 0), "positive")
  expect_error(matched_null_test(x, function(d) 1, R = 5, parallel = NA), "TRUE or FALSE")
})

test_that("every invalid return value from cluster_fn stops with an informative error", {
  set.seed(5)
  x <- matrix(rnorm(40), 20, 2)
  # each of these fails on the first call, which is the real data
  expect_error(matched_null_test(x, function(d) NA_real_,  R = 5), "returned NA on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) NaN,       R = 5), "returned NaN on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) Inf,       R = 5), "returned Inf on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) -Inf,      R = 5), "returned -Inf on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) "2",       R = 5), "class <character> and length 1 on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) factor(2), R = 5), "class <factor> and length 1 on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) TRUE,      R = 5), "class <logical> and length 1 on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) c(2, 3),   R = 5), "numeric vector of length 2 on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) NULL,      R = 5), "class <NULL> and length 0 on the real data.*finite number")
  expect_error(matched_null_test(x, function(d) list(2),   R = 5), "class <list> and length 1 on the real data.*finite number")
})

test_that("a valid real value but an invalid twin stops and names that twin", {
  set.seed(6)
  x <- matrix(rnorm(60), 30, 2)
  calls <- 0L
  flaky <- function(d) {
    calls <<- calls + 1L
    if (calls == 4L) NA_real_ else 2   # call 1 is the real data, call 4 is twin 3
  }
  expect_error(matched_null_test(x, flaky, R = 10), "returned NA on null twin 3 of 10")
  calls <- 0L
  expect_error(matched_null_test(x, flaky, R = 10), "finite number")
})

test_that("an error raised inside cluster_fn propagates unchanged", {
  set.seed(6)
  x <- matrix(rnorm(60), 30, 2)
  calls <- 0L
  boom <- function(d) {
    calls <<- calls + 1L
    if (calls == 3L) stop("pipeline failed on purpose")
    2
  }
  expect_error(matched_null_test(x, boom, R = 10), "^pipeline failed on purpose$",
               class = "simpleError")
})

test_that("an invalid twin is caught under parallel = TRUE as well", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  set.seed(6)
  x <- matrix(rnorm(60), 30, 2)
  calls <- 0L
  flaky <- function(d) {
    calls <<- calls + 1L
    if (calls == 4L) NA_real_ else 2
  }
  # future.apply also warns that it is cancelling the remaining iterations;
  # the error is what this test is about
  expect_error(suppressWarnings(matched_null_test(x, flaky, R = 10, parallel = TRUE)),
               "returned NA on null twin 3 of 10.*finite number")
})

test_that("attributes on a valid return value are stripped", {
  set.seed(9)
  x <- matrix(rnorm(40), 20, 2)
  named_int <- function(d) c(k = 3L)
  r <- matched_null_test(x, named_int, R = 5)
  expect_identical(r$real, 3)
  expect_identical(r$null, rep(3, 5))
})

# A pipeline returning a continuous summary of the twin itself, so that two runs
# agree only when they drew the same twins. A count statistic would be constant
# on data this small and the comparisons below would pass vacuously.
twin_summary <- function(d) sum(d[1:5, 1])

test_that("a seeded run reproduces exactly and reports the requested R", {
  set.seed(3)
  x <- matrix(rnorm(60 * 3), 60, 3)
  R <- 25
  set.seed(11); a <- matched_null_test(x, twin_summary, R = R)
  set.seed(11); b <- matched_null_test(x, twin_summary, R = R)
  expect_identical(a$null, b$null)
  expect_identical(a$p_exceed, b$p_exceed)
  expect_identical(a$interval, b$interval)
  expect_length(a$null, R)
  expect_identical(a$R, R)
  expect_identical(a$p_exceed, (1 + sum(a$null >= a$real)) / (R + 1))
  expect_gt(stats::var(a$null), 0)
})

test_that("parallel = TRUE is reproducible from a seed", {
  skip_if_not_installed("future.apply")
  set.seed(3)
  x <- matrix(rnorm(60 * 3), 60, 3)

  set.seed(11); a <- matched_null_test(x, twin_summary, R = 9, parallel = TRUE)$null
  set.seed(11); b <- matched_null_test(x, twin_summary, R = 9, parallel = TRUE)$null

  expect_equal(a, b)
  expect_gt(stats::var(a), 0)          # the statistic really does vary
})

test_that("parallel results do not depend on the number of workers", {
  skip_on_cran()
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  set.seed(3)
  x <- matrix(rnorm(60 * 3), 60, 3)

  future::plan(future::sequential)
  set.seed(11); one <- matched_null_test(x, twin_summary, R = 9, parallel = TRUE)$null

  future::plan(future::multisession, workers = 2)
  set.seed(11); two <- matched_null_test(x, twin_summary, R = 9, parallel = TRUE)$null

  expect_equal(one, two)
})

test_that("the printed outcome names the three cases", {
  set.seed(8)
  x <- matrix(rnorm(40), 20, 2)

  # real statistic inside a constant null
  inside <- matched_null_test(x, function(d) 2, R = 9)
  expect_true(inside$within)
  expect_output(print(inside),
                "outcome:                      no detected excess (within the reference interval)",
                fixed = TRUE)

  # real statistic above a constant null: the first call is the real data
  calls <- 0L
  above_fn <- function(d) { calls <<- calls + 1L; if (calls == 1L) 5 else 2 }
  above <- matched_null_test(x, above_fn, R = 9)
  expect_false(above$within)
  expect_gt(above$real, above$interval[2])
  expect_output(print(above),
                "departure from the reference model (above the interval; kind not settled by this test)",
                fixed = TRUE)

  # real statistic below a constant null
  calls <- 0L
  below_fn <- function(d) { calls <<- calls + 1L; if (calls == 1L) 1 else 2 }
  below <- matched_null_test(x, below_fn, R = 9)
  expect_false(below$within)
  expect_lt(below$real, below$interval[1])
  expect_output(print(below),
                "below the reference range (a diagnostic, not a verdict)",
                fixed = TRUE)
})

test_that("the printed exceedance line names the direction it counts", {
  set.seed(19)
  x <- matrix(rnorm(40), 20, 2)

  # deterministic pipelines: the first call is always the real data, every
  # twin after it answers 2. With real = 3 no twin reaches the real value, so
  # the numerator counts none of them and the exceedance is 1 / (R + 1).
  calls <- 0L
  above_fn <- function(d) { calls <<- calls + 1L; if (calls == 1L) 3 else 2 }
  above <- matched_null_test(x, above_fn, R = 9)

  # with real = 1 every twin reaches it, so the numerator counts all nine.
  calls <- 0L
  below_fn <- function(d) { calls <<- calls + 1L; if (calls == 1L) 1 else 2 }
  below <- matched_null_test(x, below_fn, R = 9)

  for (fit in list(above, below)) {
    expect_equal(fit$p_exceed, (1 + sum(fit$null >= fit$real)) / (fit$R + 1))
    expect_output(print(fit), "MC exceedance (null >= real): ", fixed = TRUE)

    line <- grep("MC exceedance", capture.output(print(fit)), value = TRUE)
    expect_length(line, 1L)
    expect_equal(as.numeric(sub(".*: +", "", line)),
                 round((1 + sum(fit$null >= fit$real)) / (fit$R + 1), 3))
  }

  # the two directions are not interchangeable, so a label naming the wrong
  # one would be reporting the complement on one of these runs
  expect_equal(above$p_exceed, 0.1)
  expect_equal(below$p_exceed, 1)
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

test_that("the seeded parallel = FALSE stream matches the 0.2.1 stream", {
  # Values recorded from matchednull 0.2.1 with this exact pipeline and seed.
  # They pin the default stream so that a later change cannot silently move it.
  set.seed(3)
  x <- matrix(rnorm(60 * 3), 60, 3)
  set.seed(11); a <- matched_null_test(x, twin_summary, R = 25)
  expect_equal(a$real, -1.95201998218585)
  expect_equal(unname(a$interval), c(-3.3755369745907, 3.46709137829001))
  expect_equal(a$p_exceed, 23 / 26)
})
