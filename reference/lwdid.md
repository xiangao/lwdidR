# Lee-Wooldridge DiD via unit-specific pre-treatment transformations

Main entry point implementing the Lee & Wooldridge (2026) panel
difference-in-differences estimator. Supports both common-timing and
staggered adoption designs. Each unit's outcome is residualised using
only its own pre-treatment observations before running a pooled
cross-sectional OLS.

## Usage

``` r
lwdid(
  data,
  y,
  ivar,
  tvar,
  gvar = NULL,
  post = NULL,
  dvar = NULL,
  rolling = "demean",
  method = NULL,
  control_group = "never_treated",
  aggregate = "overall",
  vce = NULL,
  cluster_var = NULL,
  controls = NULL,
  season_var = NULL,
  pre = -1L,
  never = FALSE,
  attgt = FALSE,
  ydot = FALSE,
  reps = 999L,
  seed = NULL,
  level = 95,
  nose = FALSE,
  nboot = 999,
  nperm = 999,
  vce_inner = "hc3"
)
```

## Arguments

- data:

  A long-format panel data frame (one row per unit-period).

- y:

  Character. Name of the outcome column.

- ivar:

  Character. Name of the unit identifier column.

- tvar:

  Character. Name of the calendar time column (numeric or integer).

- gvar:

  Character or NULL. Name of the first-treatment-year column for
  staggered designs. Units with value `0` or `NA` are treated as
  never-treated. Set to `NULL` (default) for common-timing designs
  (supply `post` instead).

- post:

  Character or NULL. Name of a binary post-treatment indicator column (0
  = pre, 1 = post). Required when `gvar = NULL`; ignored otherwise.

- dvar:

  Character or NULL. Name of a unit-level treatment-group indicator for
  common-timing designs. If omitted, `lwdid()` keeps
  backwards-compatible behavior when `post` is the treatment-on
  indicator; if `post` is a calendar post indicator, it looks for an
  unambiguous unit-invariant treatment column named `treat`, `treated`,
  `D`, or `d`.

- rolling:

  Character. Transformation method applied to each unit's pre-treatment
  observations:

  `"demean"`

  :   Subtract the unit's pre-period mean (default).

  `"detrend"`

  :   Remove a linear trend fitted on pre-periods.

  `"demeanq"`

  :   Seasonal demeaning; requires `season_var`.

  `"detrendq"`

  :   Seasonal detrending; requires `season_var`.

- method:

  Character or NULL. If one of `"ra"`, `"ipw"`, `"ipwra"`, the large-N
  path (Lee & Wooldridge 2026a) is used: per-cohort/period ATT(g,t) via
  regression adjustment, inverse-probability weighting, or doubly-robust
  IPWRA, aggregated to the event-study WATT(r) path with
  wild-cluster-bootstrap inference. `ipw`/`ipwra` require `controls`.
  `rolling` must be `"demean"` or `"detrend"`. If `NULL` (default), the
  small-N path is used.

- control_group:

  Character. Control group for staggered designs: `"never_treated"`
  (default) or `"not_yet_treated"`.

- aggregate:

  Character. Aggregation level for staggered designs: `"overall"`
  (default), `"cohort"`, or `"none"` (returns all (g,r) pairs).

- vce:

  Character or NULL. Variance-covariance estimator: `NULL`
  (homoskedastic OLS), `"hc1"`, `"hc3"`, `"cluster"`, `"wildboot"` (wild
  cluster bootstrap), or `"permutation"` (randomisation inference). The
  last two are distribution-free and recommended at small N.

- cluster_var:

  Character or NULL. Column name for clustering; required when
  `vce = "cluster"`.

- controls:

  Character vector or NULL. Names of time-invariant control variables to
  include in the cross-sectional regression.

- season_var:

  Character or NULL. Column name of the seasonal indicator (required for
  `rolling = "demeanq"` or `"detrendq"`).

- pre:

  Integer. Large-N only. Number of most-recent pre-periods used in the
  transformation; `-1` (default) uses all pre-periods.

- never:

  Logical. Large-N only. If `TRUE`, the comparison group is
  never-treated units only (default `FALSE` uses not-yet-treated units
  too).

- attgt:

  Logical. Large-N only. If `TRUE`, the ATT(g,t) cell estimates are
  returned (and printed).

- ydot:

  Logical. Large-N only. If `TRUE`, the per-cohort transformed outcomes
  are returned.

- reps:

  Integer. Large-N only. Wild cluster bootstrap replications (default
  999).

- seed:

  Integer or NULL. Large-N only. Seed for the wild bootstrap.

- level:

  Numeric. Confidence level for the large-N path (default 95).

- nose:

  Logical. Large-N only. If `TRUE`, skip standard errors (faster).

- nboot:

  Integer. Number of bootstrap replications for `vce = "wildboot"`
  (default 999).

- nperm:

  Integer. Number of permutations for `vce = "permutation"` (default
  999).

- vce_inner:

  Character. Inner variance estimator used when computing the observed
  t-statistic inside the wild bootstrap (default `"hc3"`).

## Value

An object of class `"lwdid"`, a list containing:

- `design`:

  `"staggered"` or `"common_timing"`.

- `att_overall`:

  Estimated overall ATT.

- `se_overall`:

  Standard error of overall ATT.

- `tstat`:

  t-statistic.

- `pvalue`:

  Two-sided p-value.

- `att_by_cohort`:

  Data frame of cohort-specific ATTs (staggered only).

- `att_by_cohort_time`:

  Data frame of (g,r)-specific ATTs (staggered only).

- `att_by_period`:

  Data frame of period-specific ATTs (common timing only).

- `ci_lower`, `ci_upper`:

  95% confidence interval bounds (common timing only).

- `N`:

  Sample size at first post-treatment period (common timing only).

## References

Lee, S. J., & Wooldridge, J. M. (2026a). Simple Transformation Approach
to Difference-in-Differences Estimation for Panel Data. Journal of
Business & Economic Statistics, 1-27.
[doi:10.1080/07350015.2026.2683047](https://doi.org/10.1080/07350015.2026.2683047)

Lee, S. J., & Wooldridge, J. M. (2026b). Simple Approaches to Inference
with Difference-in-Differences Estimators with Small Cross-Sectional
Sample Sizes. Working paper, SSRN 5325686.
<https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686>

`lwdidR` is an R port of the original Stata `lwdid` command by Soo Jeong
Lee and Jeffrey M. Wooldridge: <https://github.com/Soo-econ/lwdid>. The
Castle Doctrine example (Section 7.2) and the small-sample inference
options are from (2026b).

## Examples

``` r
# Load bundled Castle Doctrine dataset
castle <- read.csv(system.file("extdata", "castle.csv", package = "lwdidR"))
castle$gvar <- castle$effyear
castle$gvar[is.na(castle$gvar) | castle$gvar == 0] <- NA

# Staggered design with demeaning and HC3 standard errors
res <- lwdid(castle, "lhomicide", "sid", "year",
             gvar = "gvar", rolling = "demean", vce = "hc3")
print(res)
#> 
#> Lee-Wooldridge DiD (lwdidR)
#> Design:      staggered
#> Transf.:     demean
#> VCE:         HC3
#> --------------------------------------------------
#> Overall ATT:   0.0917
#> SE:            0.0612
#> t-stat:        1.4997
#> p-value:       0.1402
#> --------------------------------------------------
#> 
#> Cohort-specific effects:
#>  cohort     att      se tstat    pvalue
#>    2005 0.08017 0.03215 2.493 1.884e-02
#>    2006 0.06824 0.08920 0.765 4.488e-01
#>    2007 0.11406 0.09838 1.159 2.552e-01
#>    2008 0.14605 0.08203 1.780 8.548e-02
#>    2009 0.21108 0.03550 5.946 2.115e-06
#> 
```
