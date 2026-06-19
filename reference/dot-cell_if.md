# Per-cell influence function for the ATT(g,t) coefficient

Per-cell influence function for the ATT(g,t) coefficient

## Usage

``` r
.cell_if(cell, method)
```

## Arguments

- cell:

  A fitted-cell object from `.fit_cell_large` (with `$dat`, `$fit`).

- method:

  "ra", "ipw", or "ipwra".

## Value

Numeric vector of length = \#used rows of the cell (the IF on the
estimation sample).
