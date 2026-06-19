# Aggregate ATT(g,t) cells to WATT(r) + Pre/Post averages

Aggregate ATT(g,t) cells to WATT(r) + Pre/Post averages

## Usage

``` r
.aggregate_large(attgt, rolling)
```

## Arguments

- attgt:

  Data frame with columns g, t, r, att, ngt (one row per cell).

- rolling:

  "demean" or "detrend" (controls the anchor periods).

## Value

Data frame with columns effect, ryear, watt, n_cells, n_units. Rows:
Pre_avg, Post_avg, then one row per relative period r (sorted). SE / t /
p / CI columns are added later by the inference layer (Phase 2).
