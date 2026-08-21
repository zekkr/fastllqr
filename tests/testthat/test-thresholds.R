test_that("LLQR threshold bases follow the documented formulas", {
  n <- 1000
  h <- 0.2
  expect_equal(
    fastllqr:::.llqr_threshold_base(n, h, "gaussian"),
    log(log(n)) / sqrt(log(n))
  )
  expect_equal(
    fastllqr:::.llqr_threshold_base(n, h, "epanechnikov"),
    sqrt(log(n) / (n * h))
  )
})

test_that("TVCQR threshold base uses the supplied predictor norm", {
  x <- rbind(c(3, 4), c(0, 2))
  h <- 0.25
  expect_equal(
    fastllqr:::.tvcqr_threshold_base(x, h),
    log(nrow(x))^4 * h^2 * 5
  )
})

test_that("initial threshold multiplies by log log n", {
  expect_equal(
    fastllqr:::.initial_threshold(0.1, 2, 100),
    0.2 * log(log(100))
  )
})
