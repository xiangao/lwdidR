# Apply unit-specific transformation to panel data

Computes ydot (residualised outcome), ydot_postavg (post-treatment
average of ydot per unit), and marks the firstpost cross-section.

## Usage

``` r
apply_transform(df, y, ivar, tindex, post, rolling, tpost1, season_var = NULL)
```

## Arguments

- df:

  Data frame (long format, one row per unit-period)

- y:

  Outcome column name

- ivar:

  Unit identifier column name

- tindex:

  Integer time index column name

- post:

  Binary post-treatment indicator column name (0=pre, 1=post)

- rolling:

  Transformation method: "demean", "detrend", "demeanq", "detrendq"

- tpost1:

  First post-treatment period index

- season_var:

  Column name of seasonal indicator (required for demeanq/detrendq)

## Value

Data frame with added columns: ydot, ydot_postavg, firstpost
