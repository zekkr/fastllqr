.validate_numeric_vector <- function(x, name, min_length = 1L) {
  if (!is.numeric(x) || is.object(x)) {
    stop(name, " must be a plain numeric vector.", call. = FALSE)
  }
  x <- as.numeric(x)
  if (length(x) < min_length || anyNA(x) || any(!is.finite(x))) {
    stop(name, " must contain at least ", min_length,
         " finite, non-missing value(s).", call. = FALSE)
  }
  x
}

.validate_tau <- function(tau) {
  tau <- as.numeric(tau)
  if (length(tau) != 1L || is.na(tau) || !is.finite(tau) ||
      tau <= 0 || tau >= 1) {
    stop("tau must be a finite scalar strictly between 0 and 1.", call. = FALSE)
  }
  tau
}

.validate_positive_scalar <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(name, " must be a positive finite scalar.", call. = FALSE)
  }
  x
}

.validate_threshold_factor <- function(threshold_factor) {
  .validate_positive_scalar(threshold_factor, "threshold_factor")
}

.validate_tvcqr_x <- function(x) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }
  if (is.vector(x) && is.numeric(x)) {
    x <- matrix(x, ncol = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L ||
      anyNA(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric matrix with at least one row and one column.",
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}
