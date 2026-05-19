# lwdidR

[![pkgdown](https://img.shields.io/badge/pkgdown-site-blue.svg)](https://xiangao.github.io/lwdidR/)

R implementation of the Lee & Wooldridge (2025) panel
difference-in-differences estimator via unit-specific pre-treatment
transformations.

## Installation

``` r

# Install from GitHub
devtools::install_github("xiangao/lwdidR")
```

## Overview

`lwdidR` implements the estimator from:

> Lee, Y., & Wooldridge, J. M. (2025). A simple panel data approach to
> difference-in-differences under general treatment patterns.

The key idea: residualise each unit’s outcome using only its own
pre-treatment observations (demean or detrend), then run a pooled
cross-sectional OLS. This sidesteps the Callaway-Sant’Anna / Sun-Abraham
aggregation issues while allowing flexible control for unit-specific
trends.

**Supports:** - Common-timing and staggered adoption designs -
Transformation methods: `demean`, `detrend`, `demeanq`, `detrendq` -
Standard errors: homoskedastic OLS, HC1, HC3, cluster-robust -
Cohort-specific ATTs and (g,r)-level event-study estimates

## Quick Start

``` r

library(lwdidR)

# Load bundled Castle Doctrine dataset (50 US states, 2000-2010)
castle <- read.csv(system.file("extdata", "castle.csv", package = "lwdidR"))
castle$gvar <- castle$effyear
castle$gvar[is.na(castle$gvar) | castle$gvar == 0] <- NA

# Staggered design: demeaning with HC3 SEs
res <- lwdid(castle, "lhomicide", "sid", "year",
             gvar = "gvar", rolling = "demean", vce = "hc3")
print(res)
```

    Lee-Wooldridge DiD (lwdidR)
    Design:      staggered
    Transf.:     demean
    VCE:         HC3
    --------------------------------------------------
    Overall ATT:   0.0917
    SE:            0.0612
    t-stat:        1.4995
    p-value:       0.1376

## Replication: Lee & Wooldridge (2025), Section 7.2

Castle Doctrine laws and log homicide rates:

| Method              | ATT    | SE     | t-stat |
|---------------------|--------|--------|--------|
| Demeaning (OLS SE)  | 0.0917 | 0.0571 | 1.607  |
| Demeaning (HC3 SE)  | 0.0917 | 0.0612 | 1.500  |
| Detrending (HC3 SE) | 0.0666 | 0.0550 | 1.210  |

All results match paper Table 7.2 within tolerance 0.001. See
[`vignette("castle_law")`](https://xiangao.github.io/lwdidR/articles/castle_law.md)
for full replication.

## Key Functions

| Function | Description |
|----|----|
| [`lwdid()`](https://xiangao.github.io/lwdidR/reference/lwdid.md) | Main estimator (common-timing and staggered) |
| [`print.lwdid()`](https://xiangao.github.io/lwdidR/reference/print.lwdid.md) | Compact results display |
| [`summary.lwdid()`](https://xiangao.github.io/lwdidR/reference/summary.lwdid.md) | Results + period/cohort details |

## Documentation & vignettes

Full documentation: **<https://xiangao.github.io/lwdidR/>**

| Page | Description |
|----|----|
| [Castle Law replication](https://xiangao.github.io/lwdidR/articles/castle_law.html) | Full Lee-Wooldridge Section 7.2 replication |
| [Simulation comparison](https://xiangao.github.io/lwdidR/articles/simulation_comparison.html) | Simulation comparison across transformations |
| [`lwdid()`](https://xiangao.github.io/lwdidR/reference/lwdid.html) | Main estimator |
| [Reference index](https://xiangao.github.io/lwdidR/reference/index.html) | All documented functions on one page |

## Algorithm Notes

**Overall ATT (staggered):** Not a delta-method average of cohort ATTs.
Instead, a pooled cross-sectional regression: - Treated unit i: outcome
= ydot_postavg for its own cohort g - Never-treated unit i: outcome =
weighted average of ydot_postavg across all cohorts - Single OLS on
N_treated + N_never_treated units

This matches the Python reference implementation and yields SE ≈ 0.061
(not 0.051 from naive averaging).

## Reference

Lee, Y., & Wooldridge, J. M. (2025). A simple panel data approach to
difference-in-differences under general treatment patterns. *Working
paper.*
