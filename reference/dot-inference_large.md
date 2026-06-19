# Add wild-cluster-bootstrap SEs + sup-t bands to a large-N result

Add wild-cluster-bootstrap SEs + sup-t bands to a large-N result

## Usage

``` r
.inference_large(
  res,
  d,
  ivar,
  cluster_var = NULL,
  reps = 999L,
  seed = NULL,
  level = NULL
)
```

## Arguments

- res:

  Output of `lwdid_large`.

- d:

  The analysis-sample data frame (touse subset) used inside
  `lwdid_large`.

- ivar:

  Unit id column name (default cluster).

- cluster_var:

  Optional clustering column name (defaults to `ivar`).

- reps:

  Bootstrap replications (default 999).

- seed:

  Optional integer seed.

- level:

  Confidence level (default from res\$meta).

## Value

`res` with the `watt` data frame extended by se, se_analytic, t,
p_value, lower_ci, upper_ci.
