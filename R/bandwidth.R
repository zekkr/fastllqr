.llqr_bandwidth <- function(x, y, tau, h, kernel) {
  if (!is.null(h)) {
    return(.validate_positive_scalar(h, "h"))
  }

  n <- length(y)
  if (identical(kernel, "epanechnikov")) {
    return(n^(-0.2))
  }

  trim <- floor(0.2 * n)
  keep <- order(y)[seq.int(max(1L, trim), min(n, n - trim))]
  h_value <- suppressWarnings(KernSmooth::dpill(x[keep], y[keep]))
  h_value <- 1.25 * h_value *
    (tau * (1 - tau) / stats::dnorm(stats::qnorm(tau))^2)^0.2

  if (!is.finite(h_value) || h_value <= 0) {
    fallback_scale <- min(2, stats::sd(y))
    if (!is.finite(fallback_scale) || fallback_scale <= 0) {
      fallback_scale <- 1
    }
    h_value <- 1.25 * max(n^(-0.2), fallback_scale * n^(-0.2))
  }
  as.numeric(h_value)
}

.tvcqr_bandwidth <- function(n, h) {
  if (is.null(h)) {
    return(n^(-0.2))
  }
  .validate_positive_scalar(h, "h")
}
