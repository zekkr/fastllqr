.direct_llqr <- function(x, y, z, h, tau, kernel) {
  design <- cbind(1, x)
  vapply(z, function(point) {
    weights <- fastllqr:::.llqr_kernel_weights((point - x) / h, kernel)
    active <- weights > 0
    fit <- quantreg::rq.wfit(
      x = design[active, , drop = FALSE],
      y = y[active],
      weights = weights[active],
      tau = tau,
      method = "br"
    )
    fit$coefficients[1] + point * fit$coefficients[2]
  }, numeric(1))
}

test_that("LLQR matches direct weighted quantile regression", {
  skip_if_not_installed("quantreg")
  set.seed(101)
  x <- sort(stats::runif(100))
  y <- 1 + 2 * x^2 + stats::rnorm(100)
  z <- seq(0.1, 0.9, length.out = 25)

  for (kernel in c("gaussian", "epanechnikov")) {
    fit <- llqr_seq_ppro(x, y, tau = 0.5, z = z, kernel = kernel)
    reference <- .direct_llqr(x, y, z, fit$h, 0.5, kernel)
    expect_equal(fit$ll_est, reference, tolerance = 1e-8)
    expect_named(fit, c("ll_est", "d_ll_est", "z", "h"),
                 ignore.order = FALSE)
  }
})

test_that("LLQR restores unsorted evaluation-point order", {
  set.seed(102)
  x <- sort(stats::runif(80))
  y <- x + stats::rnorm(80)
  z <- c(0.8, 0.2, 0.5)
  fit <- llqr_seq_ppro(x, y, z = z, kernel = "epanechnikov")
  sorted <- llqr_seq_ppro(x, y, z = sort(z), kernel = "epanechnikov")
  expect_equal(fit$ll_est[order(z)], sorted$ll_est, tolerance = 1e-12)
  expect_identical(fit$z, z)
})

test_that("custom Epanechnikov bandwidth remains certified", {
  set.seed(103)
  x <- sort(stats::runif(100))
  y <- 1 + x + stats::rnorm(100)
  fit <- llqr_seq_ppro(x, y, z = seq(0.2, 0.8, length.out = 12),
                       h = 0.3, kernel = "epanechnikov")
  reference <- .direct_llqr(
    x, y, fit$z, fit$h, 0.5, "epanechnikov"
  )
  expect_true(all(is.finite(fit$ll_est)))
  expect_equal(fit$ll_est, reference, tolerance = 1e-8)
  expect_equal(
    fastllqr:::.llqr_threshold_base(length(x), 0.3, "epanechnikov"),
    sqrt(log(length(x)) / (length(x) * 0.3))
  )
})

test_that("LLQR threshold factors retain the certified full solution", {
  set.seed(104)
  x <- sort(stats::runif(100))
  y <- 1 + x + stats::rnorm(100)
  z <- seq(0.2, 0.8, length.out = 15)
  fits <- lapply(c(0.05, 0.1, 0.2), function(factor) {
    llqr_seq_ppro(x, y, z = z, kernel = "epanechnikov",
                  threshold_factor = factor)
  })
  expect_equal(fits[[1L]]$ll_est, fits[[2L]]$ll_est, tolerance = 1e-8)
  expect_equal(fits[[2L]]$ll_est, fits[[3L]]$ll_est, tolerance = 1e-8)
})
