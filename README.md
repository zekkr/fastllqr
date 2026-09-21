# fastllqr

Fortran-backed screened local linear and time-varying coefficient quantile
regression. Version 0.2.0 uses the shared U11 weighted-quantile-regression core.

```sh
R CMD INSTALL fastllqr_0.2.0.tar.gz
```

After publication of the fixed GitHub release:

```r
remotes::install_github("zekkr/fastllqr@v0.2.0")
```

R >= 4.1, KernSmooth and a compatible Fortran compiler are required. quantreg is
used only by validation and the separate reproduction code. See
https://github.com/zekkr/fast_llqr/tree/main/reproduction for reproduction modes,
archived sources, compiler profiles and the historical environment audit.

```r
set.seed(2026)
x <- runif(100)
y <- 1 + 2*x^2 + rnorm(100)
fit <- fastllqr::llqr_seq_ppro(x,y,kernel="epanechnikov")
attr(fit,"solver_info")
checked <- fastllqr::llqr_seq_ppro(x,y,kernel="epanechnikov",
                                 audit=TRUE,diagnostics=TRUE)
checked$diagnostics$H_seq
```

`llqr_seq_ppro()` and `tvcqr_seq_ppro()` retain their original default parameters
and four default return fields. `audit=TRUE` adds full R-side rank/interpolation
validation; it is not needed to enable the core's mandatory strict omitted-sign
checks. `diagnostics=TRUE` appends interpolation indices and 18 per-point counters.
The lightweight `solver_info` attribute is always present. Full-active recovery
uses U11 itself; no seq fallback is called. Kernel failure produces an explicit
R error with the evaluation index. Rank-deficient positive support is rejected,
including when a transported basis contains zero-weight padding rows.

U11 source fingerprints are in `inst/U11-PROVENANCE.json`. Canonical sources live
in the research repository; `reproduction/check_source_sync.py` checks the
package copy. The paper's original measurements used archived research calls,
not a retrospectively named package version. Public wrapper overhead and portable
compiler settings can make package timings differ from the archived timings.

For the local macOS/Homebrew compiler workaround and precise validation coverage,
see the reproduction documentation. The GitHub repositories were private at the
2026-09-21 check; reader access and public release are separate pending steps.
