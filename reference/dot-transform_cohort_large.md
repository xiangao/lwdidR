# Large-N rolling transformation for one cohort

Large-N rolling transformation for one cohort

## Usage

``` r
.transform_cohort_large(df, y, ivar, tvar, g, rolling, pre = -1L)
```

## Arguments

- df:

  Long panel data frame (one row per unit-period).

- y:

  Outcome column name.

- ivar:

  Unit id column name.

- tvar:

  Calendar time column name (integer-valued).

- g:

  Cohort (first-treatment period) to transform for.

- rolling:

  "demean" or "detrend".

- pre:

  Integer. Number of most-recent pre-periods to use; -1 uses all
  pre-periods (t \< g). With "detrend", an effective window of \>= 2 is
  required.

## Value

Numeric vector aligned to the rows of `df`: the transformed outcome
`ydot_g` (NA where the unit has no usable pre-window), with anchor cells
set to exactly 0.
