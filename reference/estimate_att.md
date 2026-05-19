# Estimate ATT via OLS on firstpost cross-section (uses ydot_postavg)

Estimate ATT via OLS on firstpost cross-section (uses ydot_postavg)

## Usage

``` r
estimate_att(
  df,
  d,
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

  Transformed data with ydot_postavg and firstpost columns

- d:

  Treatment indicator column

- vce:

  Variance estimator

- cluster_var:

  Cluster variable name

- controls:

  Control variable names

## Value

Named list with att, se, tstat, pvalue, ci_lower, ci_upper, N, df_resid
