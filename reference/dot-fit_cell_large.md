# Fit one large-N ATT(g,t) cell

Fit one large-N ATT(g,t) cell

## Usage

``` r
.fit_cell_large(ydot, dvar, xc, xraw, method)
```

## Arguments

- ydot:

  Numeric vector: cohort-g transformed outcome on the cell
  cross-section.

- dvar:

  Numeric 0/1 treatment indicator (1 = first-treated in cohort g).

- xc:

  Matrix of cohort-centered covariates (cell rows x p), or NULL.

- xraw:

  Matrix of RAW covariates for the propensity score (cell rows x p), or
  NULL.

- method:

  "ra", "ipw", or "ipwra".

## Value

list(ok, att, ngt, fit, w, used) where `used` is the logical vector of
rows entering estimation. `ok = FALSE` marks an unusable cell (skipped,
as in Stata's `cap ... continue`).
