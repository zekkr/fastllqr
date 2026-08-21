.llqr_kernel_weights <- function(u, kernel) {
  if (identical(kernel, "gaussian")) {
    return(stats::dnorm(u))
  }
  weights <- numeric(length(u))
  inside <- abs(u) <= 1
  weights[inside] <- 0.75 * (1 - u[inside]^2)
  weights
}

.certify_H <- function(beta, H, design, y, residual_tolerance = 1e-6,
                       rank_tolerance = 1e-10) {
  p <- ncol(design)
  valid_indices <- length(H) == p && !anyNA(H) &&
    all(H >= 1L & H <= nrow(design)) && !anyDuplicated(H)
  if (!valid_indices || length(beta) != p || any(!is.finite(beta))) {
    return(FALSE)
  }
  residual <- as.numeric(y - design %*% beta)
  qr(design[H, , drop = FALSE], tol = rank_tolerance)$rank == p &&
    all(abs(residual[H]) <= residual_tolerance)
}

.certify_llqr_result <- function(result, x, y, z, h, kernel) {
  design <- cbind(1, x)
  for (index in seq_along(z)) {
    weights <- .llqr_kernel_weights((z[index] - x) / h, kernel)
    active <- weights > 0
    if (sum(active) < 2L || qr(design[active, , drop = FALSE])$rank < 2L) {
      stop("The full local LLQR problem is rank deficient at evaluation point ",
           index, ".", call. = FALSE)
    }
    beta <- c(result$ll_est[index] - z[index] * result$d_ll_est[index],
              result$d_ll_est[index])
    if (!.certify_H(beta, result$H_seq[index, ], design, y)) {
      stop("The LLQR backend returned an uncertified basis at evaluation point ",
           index, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.certify_tvcqr_result <- function(result, x, y, h) {
  n <- nrow(x)
  time_index <- seq_len(n) / n
  base_design <- cbind(1, x)
  design <- cbind(base_design, base_design * time_index)
  q <- ncol(design)

  for (index in seq_len(n)) {
    weights <- 0.75 * (1 - ((index / n - time_index) / h)^2) *
      (abs(index / n - time_index) <= h)
    active <- weights > 0
    if (sum(active) < q || qr(design[active, , drop = FALSE])$rank < q) {
      stop("The full local TVCQR problem is rank deficient at time index ",
           index, ".", call. = FALSE)
    }
    if (!.certify_H(result$beta_full_est[index, ], result$H_seq[index, ],
                    design, y)) {
      stop("The TVCQR backend returned an uncertified basis at time index ",
           index, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}
