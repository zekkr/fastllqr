# U11 core shared with the research repository; see inst/U11-PROVENANCE.json.
.u11_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  value
}

.ppro_error_message <- function(model, ierr, failed_eval) {
  reason <- switch(as.character(ierr),
    `1` = "full-sample certification failed", `2` = "simplex did not converge",
    `3` = "invalid basis or rank deficient local design",
    `4` = "singular basis or non-finite numerical state",
    `5` = "an internal invariant was violated", "unknown backend error")
  paste0(model, " U11 failed at evaluation point ", failed_eval, ": ", reason,
         " (ierr=", ierr, ").")
}

.u11_backend <- function(design, y, grid, coordinate, kernel, h, tau, threshold) {
  n <- length(y); q <- ncol(design); m <- length(grid)
  if (!is.finite(threshold) || threshold <= 0)
    stop("The initial threshold must be positive and finite.", call. = FALSE)
  out <- .Fortran(C_ssqr_kernel_path,
    a=as.double(design), y=as.double(y), n=as.integer(n), q=as.integer(q),
    grid=as.double(grid), ne=as.integer(m), coordinate=as.double(coordinate),
    kernel=as.integer(kernel), h=as.double(h), tau=as.double(tau),
    tol=1e-14, maxit=1000000L, threshold=rep(as.double(threshold),m),
    min_keep=1L, cache_flags=27L, provider_flags=1L,
    beta=double(m*q), hseq=integer(m*q), diagnostics=integer(m*18L),
    ierr=0L, failed_eval=0L)
  if (out$ierr != 0L) stop(.ppro_error_message("Screened solver", out$ierr,
                                            out$failed_eval), call. = FALSE)
  beta <- matrix(out$beta,m,q); H <- matrix(out$hseq,m,q)
  if (any(!is.finite(beta)) || any(H < 1L | H > n) ||
      any(apply(H,1L,anyDuplicated) != 0L))
    stop("U11 returned an invalid solution or basis.", call. = FALSE)
  # Positive-weight basis rows certify full support rank. Zero-weight padding
  # does not: audit only those windows, without disabling the core sign checks.
  basis_coordinates <- matrix(coordinate[H],m,q)
  u <- (grid-basis_coordinates)/h
  basis_weights <- if(kernel==1L) stats::dnorm(u) else pmax(0.75*(1-u*u),0)
  needs_support_check <- which(rowSums(basis_weights > 0) < q)
  for(j in needs_support_check) {
    u <- (grid[j]-coordinate)/h
    weights <- if(kernel==1L) stats::dnorm(u) else pmax(0.75*(1-u*u),0)
    active <- weights > 0
    if(sum(active)<q || qr(design[active,,drop=FALSE])$rank < q)
      stop("The full local problem is rank deficient at evaluation point ",j,".",call.=FALSE)
  }
  d <- matrix(out$diagnostics,m,18L,dimnames=list(NULL,c(
    "n_active","first_tableau_rows","final_tableau_rows","iterations",
    "repairs","init_mode","init_trigger","full_recovery","independent",
    "certificate_hits","residual_rows","first_full_m_recovery",
    "initial_threshold_hits","threshold_expansion_steps","effective_threshold_hits",
    "basis_forced_rows","first_aggregate_rows","first_retained_size")))
  d[1L,13:18] <- NA_integer_
  info <- list(backend="unified_u11", version="0.2.0", ierr=0L, failed_eval=0L,
    fallback_triggered=FALSE, fallback_reason=NULL, returned_backend="unified_u11",
    full_active_recovery_count=sum(d[,"full_recovery"]),
    independent_init_count=sum(d[,"independent"]), repair_count=sum(d[,"repairs"]),
    threshold_initial=threshold, cache_flags=27L, provider_flags=1L)
  list(beta=beta,H_seq=H,diagnostics=d,solver_info=info,
       full_active_recovery_count=info$full_active_recovery_count)
}

.llqr_backend <- function(x, y, tau, z, h, kernel, threshold_factor) {
  threshold <- .initial_threshold(threshold_factor,.llqr_threshold_base(length(y),h,kernel),length(y))
  out <- .u11_backend(cbind(1,x), y,z,x,if(kernel=="gaussian") 1L else 2L,h,tau,threshold)
  out$ll_est <- out$beta[,1L]+z*out$beta[,2L]
  out$d_ll_est <- out$beta[,2L]
  out
}

.tvcqr_backend <- function(x, y, tau, h, threshold_factor) {
  n <- nrow(x); time <- seq_len(n)/n; X <- cbind(1,x); p <- ncol(X)
  threshold <- .initial_threshold(threshold_factor,.tvcqr_threshold_base(x,h),n)
  out <- .u11_backend(cbind(X,X*time),y,time,time,2L,h,tau,threshold)
  out$theta_ll_est <- out$beta[,seq_len(p),drop=FALSE]+out$beta[,p+seq_len(p),drop=FALSE]*time
  out$beta_full_est <- out$beta
  out
}

.u11_result <- function(result, backend, diagnostics, order=NULL, audit=FALSE) {
  info <- backend$solver_info; info$audit <- audit
  attr(result,"solver_info") <- info
  if (diagnostics) {
    if (is.null(order)) order <- seq_len(nrow(backend$H_seq))
    result$diagnostics <- list(H_seq=backend$H_seq[order,,drop=FALSE],
      path=backend$diagnostics[order,,drop=FALSE],solver_info=info)
  }
  result
}
