# Changelog

## lwdidR 0.1.0

Initial release.

- Implements Lee & Wooldridge (2025) DiD estimator via unit-specific
  pre-treatment transformations
- Supports common-timing and staggered adoption designs
- Transformation methods: `demean`, `detrend`, `demeanq`, `detrendq`
- Standard errors: homoskedastic OLS, HC1, HC3, cluster-robust (all
  manual sandwich)
- Staggered design: overall pooled ATT, cohort-specific ATT, and
  (g,r)-level effects
- Bundles Castle Doctrine dataset for replication of paper Section 7.2
- Vignette replicating all numerical results from Lee &
  Wooldridge (2025) Table 7.2
