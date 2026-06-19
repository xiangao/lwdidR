# Large-N LWDID (point estimates)

Large-N LWDID (point estimates)

## Usage

``` r
lwdid_large(
  data,
  y,
  ivar,
  tvar,
  gvar,
  rolling = "demean",
  method = "ra",
  controls = NULL,
  pre = -1L,
  never = FALSE,
  level = 95
)
```

## Arguments

- data:

  Long-format panel data frame.

- y:

  Outcome column name.

- ivar:

  Unit id column name.

- tvar:

  Calendar time column name (integer-valued).

- gvar:

  First-treatment-period column (0 or NA = never treated).

- rolling:

  "demean" or "detrend".

- method:

  "ra", "ipw", or "ipwra".

- controls:

  Character vector of covariate names, or NULL.

- pre:

  Integer; number of most-recent pre-periods used in the transform (-1 =
  all pre-periods).

- never:

  Logical; if TRUE, controls are never-treated units only.

- level:

  Confidence level (default 95).

## Value

list(attgt = cell-level data frame, watt = aggregated data frame, cells
= list of per-cell fit objects for Phase-2 inference, meta = list).
