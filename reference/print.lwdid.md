# Print an lwdid object

Displays a compact summary of estimation results including the design
type, transformation method, VCE, overall ATT, SE, t-statistic, and
p-value. For staggered designs also prints cohort-specific effects.

## Usage

``` r
# S3 method for class 'lwdid'
print(x, ...)
```

## Arguments

- x:

  An object of class `"lwdid"` as returned by
  [`lwdid()`](https://xiangao.github.io/lwdidR/reference/lwdid.md).

- ...:

  Currently ignored.

## Value

Invisibly returns `x`.
