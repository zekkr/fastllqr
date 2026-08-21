.ppro_error_message <- function(model, ierr, failed_eval) {
  reason <- switch(
    as.character(ierr),
    `1` = "full-sample certification failed",
    `2` = "simplex did not converge",
    `3` = "the interpolation basis was invalid",
    `4` = "the interpolation basis was singular",
    `5` = "an internal invariant was violated",
    "an unknown backend error occurred"
  )
  paste0(model, " ppro failed at evaluation point ", failed_eval,
         ": ", reason, " (ierr=", ierr, ").")
}

.llqr_backend <- function(x, y, tau, z, h, kernel, threshold_factor) {
  n <- length(y)
  rounds <- length(z)
  case <- if (identical(kernel, "gaussian")) 1L else 2L

  result <- .Fortran(
    C_llqr_ppro_fortran,
    x = as.double(x),
    y = as.double(y),
    z = as.double(z),
    m = as.integer(n),
    nvar = 1L,
    rounds = as.integer(rounds),
    tau = as.double(tau),
    h = as.double(h),
    tol = 1e-14,
    maxit = 1000000L,
    Mm_factor = as.double(threshold_factor),
    case_int = case,
    bland_int = 0L,
    min_subsample_size_in = 1L,
    ll_est = double(rounds),
    d_ll_est = double(rounds),
    H_mat = integer(rounds * 2L),
    full_active_recovery_count = integer(1L),
    ierr = integer(1L),
    failed_eval = integer(1L),
    always_same_h_refit_int = 1L,
    threshold_lower_bound_int = 0L,
    threshold_scale_mode_int = 1L
  )

  ierr <- as.integer(result$ierr)
  if (ierr != 0L) {
    stop(.ppro_error_message("LLQR", ierr, as.integer(result$failed_eval)),
         call. = FALSE)
  }

  list(
    ll_est = as.numeric(result$ll_est),
    d_ll_est = as.numeric(result$d_ll_est),
    H_seq = matrix(as.integer(result$H_mat), nrow = rounds, ncol = 2L),
    full_active_recovery_count = as.integer(result$full_active_recovery_count)
  )
}

.tvcqr_backend <- function(x, y, tau, h, threshold_factor) {
  n <- nrow(x)
  nvar <- ncol(x)
  q <- 2L * (nvar + 1L)

  result <- .Fortran(
    C_tvcqr_seq_ppro_fortran,
    x = as.double(x),
    y = as.double(y),
    m = as.integer(n),
    nvar = as.integer(nvar),
    tau = as.double(tau),
    h = as.double(h),
    h_factor = 1,
    tol = 1e-14,
    maxit = 1000000L,
    bland_int = 0L,
    Mm_factor = as.double(threshold_factor),
    theta_ll_est = double(n * (nvar + 1L)),
    beta_full_est = double(n * q),
    H_seq = integer(n * q),
    full_active_recovery_count = integer(1L),
    ierr = integer(1L),
    failed_eval = integer(1L),
    min_subsample_size_in = 1L,
    always_same_h_refit_int = 1L,
    threshold_lower_bound_int = 0L,
    threshold_scale_mode_int = 1L
  )

  ierr <- as.integer(result$ierr)
  if (ierr != 0L) {
    stop(.ppro_error_message("TVCQR", ierr, as.integer(result$failed_eval)),
         call. = FALSE)
  }

  list(
    theta_ll_est = matrix(result$theta_ll_est, nrow = n, ncol = nvar + 1L),
    beta_full_est = matrix(result$beta_full_est, nrow = n, ncol = q),
    H_seq = matrix(as.integer(result$H_seq), nrow = n, ncol = q),
    full_active_recovery_count = as.integer(result$full_active_recovery_count)
  )
}
