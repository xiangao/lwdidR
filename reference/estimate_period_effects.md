# Estimate period-specific ATTs (uses ydot)

Estimate period-specific ATTs (uses ydot)

## Usage

``` r
estimate_period_effects(
  df,
  d,
  tindex,
  post_periods,
  vce = NULL,
  cluster_var = NULL,
  controls = NULL,
  nboot = 999,
  nperm = 999,
  vce_inner = "hc3"
)
```

## Arguments

- df:

  Transformed panel data

- d:

  Treatment indicator column

- tindex:

  Time index column

- post_periods:

  Integer vector of post-period indices

- vce:

  Variance estimator

- cluster_var:

  Cluster variable name

- controls:

  Control variable names

## Value

Data frame with: tindex, att, se, tstat, pvalue, ci_lower, ci_upper, N
