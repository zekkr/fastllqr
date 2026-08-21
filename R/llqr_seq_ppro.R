#' Fast screened local linear quantile regression
#'
#' Fits one-dimensional local linear quantile regression along an ordered
#' evaluation grid using the exact screened warm-start ppro solver.
#'
#' @param x Numeric predictor vector.
#' @param y Numeric response vector with the same length as `x`.
#' @param tau Quantile level strictly between zero and one.
#' @param z Numeric evaluation points. Defaults to `x`.
#' @param h Positive bandwidth. When `NULL`, a kernel-specific default is used.
#' @param kernel Either `"gaussian"` or `"epanechnikov"`.
#' @param threshold_factor Positive finite multiplier for the initial screening
#'   threshold. The default is `0.1`.
#'
#' @return A list containing `ll_est`, `d_ll_est`, `z`, and `h`.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(80))
#' y <- 1 + 2 * x^2 + rnorm(80)
#' fit <- llqr_seq_ppro(x, y, z = seq(0.1, 0.9, length.out = 20),
#'                      kernel = "epanechnikov")
#' head(fit$ll_est)
llqr_seq_ppro <- function(x, y, tau = 0.5, z = NULL, h = NULL,
                          kernel = c("gaussian", "epanechnikov"),
                          threshold_factor = 0.1) {
  x <- .validate_numeric_vector(x, "x", min_length = 5L)
  y <- .validate_numeric_vector(y, "y", min_length = 5L)
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  tau <- .validate_tau(tau)
  kernel <- match.arg(kernel)
  threshold_factor <- .validate_threshold_factor(threshold_factor)
  if (is.null(z)) {
    z <- x
  } else {
    z <- .validate_numeric_vector(z, "z")
  }
  h <- .llqr_bandwidth(x, y, tau, h, kernel)
  sorted_index <- order(z, method = "radix")
  restore_index <- order(sorted_index)
  z_sorted <- z[sorted_index]

  backend <- .llqr_backend(x, y, tau, z_sorted, h, kernel,
                           threshold_factor)
  .certify_llqr_result(backend, x, y, z_sorted, h, kernel)

  list(
    ll_est = backend$ll_est[restore_index],
    d_ll_est = backend$d_ll_est[restore_index],
    z = z,
    h = h
  )
}
