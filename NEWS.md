# fastllqr 0.2.0

- Replace separate model-specific backends with the archived shared U11 core.
- Keep existing public defaults and default return fields; add `audit` and
  `diagnostics` options plus a `solver_info` provenance attribute.
- Support arbitrary bandwidths and restore evaluation order in public adapters.
- Reject deficient positive-weight support hidden by zero-weight basis padding.
- Preserve strict core verification and distinguish full-active recovery from
  seq fallback. No seq fallback is implemented by this backend.
- Include source fingerprints and links to fixed-source reproduction materials.

# fastllqr 0.1.0

Initial package release with separate LLQR/TVCQR backends.
