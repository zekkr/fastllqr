#' Fast screened time-varying coefficient quantile regression
#'
#' Fits time-varying coefficient quantile regression on the ordered grid
#' `seq_len(n) / n` using the exact screened warm-start ppro solver.
#'
#' @param x Finite numeric predictor matrix. Rows must be in time order.
#' @param y Numeric response vector with length `nrow(x)`.
#' @param tau Quantile level strictly between zero and one.
#' @param h Positive bandwidth. Defaults to `nrow(x)^(-1/5)`.
#' @param threshold_factor Positive finite multiplier for the initial screening
#'   threshold. The default is `1e-5`.
#'
#' @param audit Run an additional full R-side rank and interpolation audit.
#' @param diagnostics Append the interpolation indices and per-point solver diagnostics.
#' @return A list containing `theta_ll_est`, `beta_full_est`, `time_index`,
#'   and `h`.
#' The `solver_info` attribute records backend and recovery/fallback status.
#' With `diagnostics=TRUE`, an additional diagnostics element is returned.
#' @export
#'
#' @examples
#' set.seed(2)
#' n <- 80
#' time_index <- seq_len(n) / n
#' x <- matrix(rnorm(n * 2), ncol = 2)
#' y <- sin(2 * pi * time_index) + 0.5 * x[, 1] + rnorm(n)
#' fit <- tvcqr_seq_ppro(x, y)
#' head(fit$theta_ll_est)
tvcqr_seq_ppro <- function(x, y, tau = 0.5, h = NULL,
                           threshold_factor = 1e-5, audit = FALSE, diagnostics = FALSE) {
  audit <- .u11_flag(audit, "audit")
  diagnostics <- .u11_flag(diagnostics, "diagnostics")
  x <- .validate_tvcqr_x(x)
  y <- .validate_numeric_vector(y, "y")
  n <- nrow(x)
  if (length(y) != n) {
    stop("length(y) must equal nrow(x).", call. = FALSE)
  }
  q <- 2L * (ncol(x) + 1L)
  if (n < q) {
    stop("x has too few rows for the TVCQR local linear design.", call. = FALSE)
  }
  tau <- .validate_tau(tau)
  threshold_factor <- .validate_threshold_factor(threshold_factor)
  h <- .tvcqr_bandwidth(n, h)

  backend <- .tvcqr_backend(x, y, tau, h, threshold_factor)
  if (audit) .certify_tvcqr_result(backend, x, y, h)

  result <- list(
    theta_ll_est = backend$theta_ll_est,
    beta_full_est = backend$beta_full_est,
    time_index = seq_len(n) / n,
    h = h
  )
  .u11_result(result, backend, diagnostics, audit=audit)
}
