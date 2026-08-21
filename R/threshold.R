.threshold_loglog <- function(n) {
  value <- log(log(n))
  if (!is.finite(value) || value <= 0) {
    stop("The sample size is too small for the log(log(n)) threshold.",
         call. = FALSE)
  }
  value
}

.llqr_threshold_base <- function(n, h, kernel) {
  if (identical(kernel, "gaussian")) {
    return(.threshold_loglog(n) / sqrt(log(n)))
  }
  sqrt(log(n) / (n * h))
}

.tvcqr_threshold_base <- function(x, h) {
  n <- nrow(x)
  max_norm <- max(sqrt(rowSums(x^2)))
  log(n)^4 * h^2 * max_norm
}

.initial_threshold <- function(threshold_factor, threshold_base, n) {
  threshold_factor * threshold_base * .threshold_loglog(n)
}
