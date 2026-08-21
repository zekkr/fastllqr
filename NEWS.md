# fastllqr 0.1.0

- Added `llqr_seq_ppro()` for Gaussian- and Epanechnikov-kernel LLQR.
- Added `tvcqr_seq_ppro()` for ordered TVCQR fits.
- Added model-specific screening threshold profiles with public
  `threshold_factor` controls.
- Removed sequential-solver fallback and user-facing debug diagnostics.
- Added registered Fortran routines, full candidate checks, tests, and
  cross-platform R package checks.
