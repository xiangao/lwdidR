# LWDID for staggered adoption (gvar case)

Python-equivalent algorithm: For each cohort g:

1.  Apply transformation using pre-periods (t \< g) for ALL units

2.  Compute ydot_postavg(g,i) per unit Overall ATT: pooled regression

- Treated unit i: Y_bar_i = ydot_postavg(g_i, i)

- Never-treated unit i: Y_bar_i = sum_g w_g \* ydot_postavg(g,i)

- OLS: Y_bar ~ D (D=1 if ever treated) on N_treat + N_NT units Cohort
  ATT_g: per-cohort regression

- Y_bar(g,i) = ydot_postavg(g,i) for cohort g units + NT controls

- OLS: Y_bar_g ~ D_g

## Usage

``` r
lwdid_staggered(
  data,
  y,
  ivar,
  tvar,
  gvar,
  rolling = "demean",
  control_group = "never_treated",
  aggregate = "overall",
  vce = NULL,
  cluster_var = NULL,
  controls = NULL,
  season_var = NULL,
  nboot = 999,
  nperm = 999,
  vce_inner = "hc3"
)
```

## Arguments

- data:

  Long-format panel data

- y:

  Outcome column

- ivar:

  Unit identifier

- tvar:

  Calendar time column

- gvar:

  First treatment year (0/NA = never treated)

- rolling:

  Transformation method

- control_group:

  "never_treated" or "not_yet_treated"

- aggregate:

  "overall", "cohort", or "none"

- vce:

  Variance estimator

- cluster_var:

  Cluster variable

- controls:

  Control variable names

- season_var:

  Seasonal indicator column

## Value

Named list with att_overall, att_by_cohort, att_by_cohort_time
