.direct_tvcqr <- function(x, y, h, tau, indices) {
  n <- nrow(x)
  time_index <- seq_len(n) / n
  base_design <- cbind(1, x)
  design <- cbind(base_design, base_design * time_index)
  p <- ncol(base_design)

  t(vapply(indices, function(index) {
    weights <- 0.75 * (1 - ((index / n - time_index) / h)^2) *
      (abs(index / n - time_index) <= h)
    active <- weights > 0
    fit <- quantreg::rq.wfit(
      x = design[active, , drop = FALSE],
      y = y[active],
      weights = weights[active],
      tau = tau,
      method = "br"
    )
    beta <- as.numeric(fit$coefficients)
    beta[seq_len(p)] + (index / n) * beta[p + seq_len(p)]
  }, numeric(p)))
}

test_that("TVCQR matches direct weighted quantile regression", {
  skip_if_not_installed("quantreg")
  set.seed(201)
  n <- 80
  time_index <- seq_len(n) / n
  x <- matrix(stats::rnorm(n * 2), ncol = 2)
  y <- sin(2 * pi * time_index) + 0.5 * x[, 1] + stats::rnorm(n)

  fit <- tvcqr_seq_ppro(x, y)
  indices <- c(1L, 20L, 40L, 60L, 80L)
  reference <- .direct_tvcqr(x, y, fit$h, 0.5, indices)
  expect_equal(fit$theta_ll_est[indices, , drop = FALSE], reference,
               tolerance = 1e-8)
  expect_identical(fit$time_index, time_index)
  expect_named(fit,
               c("theta_ll_est", "beta_full_est", "time_index", "h"),
               ignore.order = FALSE)
})

test_that("TVCQR threshold factors do not expose diagnostics", {
  set.seed(202)
  n <- 60
  x <- matrix(stats::rnorm(n * 2), ncol = 2)
  y <- stats::rnorm(n)
  low <- tvcqr_seq_ppro(x, y, threshold_factor = 0.5e-5)
  default <- tvcqr_seq_ppro(x, y, threshold_factor = 1e-5)
  high <- tvcqr_seq_ppro(x, y, threshold_factor = 2e-5)
  expect_equal(low$theta_ll_est, default$theta_ll_est, tolerance = 1e-8)
  expect_equal(default$theta_ll_est, high$theta_ll_est, tolerance = 1e-8)
  expect_false(any(c("diagnostics", "H_seq", "returned_backend") %in% names(low)))
})
