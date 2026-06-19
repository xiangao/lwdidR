# lwdidR 0.2.0

* Added the **large-N path** of Lee & Wooldridge (2026a), activated via
  `method = "ra"`, `"ipw"`, or `"ipwra"`. It estimates ATT(g,t) by cohort and
  period (regression adjustment / inverse-probability weighting / doubly-robust
  IPWRA), aggregates to the event-study WATT(r) path and Pre/Post averages, and
  reports wild cluster bootstrap standard errors with simultaneous (sup-t)
  confidence bands. New options: `method`, `pre`, `never`, `attgt`, `ydot`,
  `reps`, `seed`, `level`, `nose`.
* Added `plot()` (event-study graph via ggplot2) and `print()`/`summary()`
  methods for the large-N result (class `lwdid_large`).
* Influence-function machinery ported verbatim from Stata `lwdid` v2.4.2; the
  large-N estimates and inference are validated against Stata on `lw_smoking` and
  `lw_walmart` (see `data-raw/stata-reference/`).

# lwdidR 0.1.2

* Corrected citation to the published version: Lee, S. J., & Wooldridge, J. M.
  (2026), "Simple Transformation Approach to Difference-in-Differences Estimation
  for Panel Data", *Journal of Business & Economic Statistics*, 1–27,
  <doi:10.1080/07350015.2026.2683047>. Fixed the first author's name (previously
  shown incorrectly as "Lee, Y.").
* Added `inst/CITATION` and acknowledged the original Stata `lwdid` package
  (<https://github.com/Soo-econ/lwdid>) in the README and function documentation.
* Corrected paper attribution for the Castle Doctrine example and the small-sample
  inference options (wild bootstrap, permutation): these come from the small-N
  paper Lee & Wooldridge (2026b), "Simple Approaches to Inference with
  Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes"
  (SSRN 5325686), Section 7.2 — not the JBES transformation paper. The package now
  cites both papers. Removed a reference to a non-existent "Table 7.2" (the Castle
  results appear in the text of Section 7.2) and corrected the cohort-effect
  equation reference to (7.10).

# lwdidR 0.1.0

Initial release.

* Implements Lee & Wooldridge (2026a, 2026b) DiD estimator via unit-specific pre-treatment transformations
* Supports common-timing and staggered adoption designs
* Transformation methods: `demean`, `detrend`, `demeanq`, `detrendq`
* Standard errors: homoskedastic OLS, HC1, HC3, cluster-robust (all manual sandwich)
* Staggered design: overall pooled ATT, cohort-specific ATT, and (g,r)-level effects
* Bundles Castle Doctrine dataset for replication of paper Section 7.2
* Vignette replicating the Castle Doctrine results in Lee & Wooldridge (2026b) Section 7.2
