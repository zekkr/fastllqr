# fastllqr

`fastllqr` provides Fortran-backed exact screened warm-start solvers for:

- one-dimensional local linear quantile regression (LLQR), and
- time-varying coefficient quantile regression (TVCQR).

The package implements the fixed algorithm profile used for the paper
*A fast simplex algorithm for local linear quantile regression*.  It does not
fall back to a sequential solver.  Every screened solution is checked against
the omitted observations; an uncertifiable reduced problem is expanded or
handled by the ppro solver's internal full-active recovery.  A terminal failure
is reported as an R error.

## Installation

```r
remotes::install_github("zekkr/fastllqr")
```

Installation from source requires a working Fortran compiler.  On Windows use
the Rtools version matching R.  On macOS or Linux install a `gfortran` version
compatible with the installed R toolchain.

## LLQR

```r
library(fastllqr)

set.seed(1)
x <- sort(runif(200))
y <- 1 + 2 * x^2 + rnorm(200)

fit <- llqr_seq_ppro(
  x,
  y,
  z = seq(0.1, 0.9, length.out = 50),
  kernel = "epanechnikov"
)

head(fit$ll_est)
```

The result contains `ll_est`, `d_ll_est`, `z`, and `h`.

## TVCQR

```r
set.seed(2)
n <- 200
time_index <- seq_len(n) / n
x <- matrix(rnorm(n * 2), ncol = 2)
y <- sin(2 * pi * time_index) + 0.5 * x[, 1] + rnorm(n)

fit <- tvcqr_seq_ppro(x, y)
head(fit$theta_ll_est)
```

The result contains `theta_ll_est`, `beta_full_est`, `time_index`, and `h`.
Rows of `x` and entries of `y` must already be in time order.

## Screening thresholds

The initial threshold has the form

\[
\gamma_n = c_\gamma\,\bar\pi_n\,\log\log n,
\]

where `threshold_factor` is \(c_\gamma\).  The fixed package profiles are:

| Model | Default `threshold_factor` | Threshold base \(\bar\pi_n\) |
|---|---:|---:|
| LLQR, Gaussian kernel | `0.1` | \(\log\log(n)/\sqrt{\log(n)}\) |
| LLQR, Epanechnikov kernel | `0.1` | \(\sqrt{\log(n)/(nh)}\) |
| TVCQR | `1e-5` | \((\log n)^4h^2\max_i\lVert x_i\rVert\) |

The Gaussian LLQR profile is an empirical full-order-grid rule.  The
Epanechnikov rule uses the actual supplied bandwidth.  The TVCQR norm is taken
over the predictor columns supplied by the user, excluding the internally
added intercept.

Increasing `threshold_factor` retains a larger initial problem; decreasing it
screens more aggressively.  Full omitted-sign verification remains mandatory
for every setting.

## Reproducibility and provenance

The initial implementation was extracted from
[`zekkr/fast_llqr`](https://github.com/zekkr/fast_llqr) commit
`98592831d5e834b4379efa0b2d84b9b876763a1b`.  Only the two ppro solvers and
their required package infrastructure are included here.  This package is
currently distributed through GitHub and is not submitted to CRAN.

## License

GPL-3.0-or-later.
