# lwdidR

[![pkgdown](https://img.shields.io/badge/pkgdown-site-blue.svg)](https://xiangao.github.io/lwdidR/)

`lwdidR` is an R implementation of the Lee and Wooldridge (2026a, 2026b)
panel difference-in-differences estimator. The estimator first
transforms outcomes using only each unit’s pre-treatment observations,
then estimates the treatment effect in a cross section.

## Installation

``` r

# Install from GitHub
devtools::install_github("xiangao/lwdidR")
```

## Overview

`lwdidR` implements the rolling difference-in-differences method of Lee
and Wooldridge, drawing on two companion papers:

> Lee, S. J., & Wooldridge, J. M. (2026a). Simple Transformation
> Approach to Difference-in-Differences Estimation for Panel Data.
> *Journal of Business & Economic Statistics*, 1–27.
> <https://doi.org/10.1080/07350015.2026.2683047>
>
> Lee, S. J., & Wooldridge, J. M. (2026b). Simple Approaches to
> Inference with Difference-in-Differences Estimators with Small
> Cross-Sectional Sample Sizes. Working paper, SSRN 5325686.
> <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686>

The outcome transformation underlying the estimator is introduced in
(2026a); the small-sample inference options (wild cluster bootstrap,
permutation) and the staggered Castle Doctrine application replicated
here come from (2026b). `lwdidR` is an R port of the original Stata
`lwdid` command by Soo Jeong Lee and Jeffrey M. Wooldridge
(<https://github.com/Soo-econ/lwdid>), whose v1.0 is based on (2026b).

The key idea is simple: residualize each unit’s outcome using its own
pre-treatment observations, either by demeaning or detrending, and then
run a pooled cross-sectional OLS.

The package supports: - Common-timing and staggered adoption designs -
Transformation methods: `demean`, `detrend`, `demeanq`, `detrendq` -
Standard errors: homoskedastic OLS, HC1, HC3, cluster-robust, wild
cluster bootstrap, permutation - Cohort-specific ATTs and (g,r)-level
event-study estimates - Large-N estimators
(`method = "ra"/"ipw"/"ipwra"`) with WATT(r) event study, wild cluster
bootstrap, and sup-t bands

### Large-N path (Lee & Wooldridge 2026a)

Setting `method` to `"ra"`, `"ipw"`, or `"ipwra"` activates the large-N
path: it estimates ATT(g,t) by cohort and period (regression adjustment,
inverse-probability weighting, or doubly-robust IPWRA), aggregates to
the event-study **WATT(r)** path plus Pre/Post averages, and reports
**wild cluster bootstrap** standard errors with **simultaneous (sup-t)
confidence bands**.

``` r

# doubly-robust IPWRA, detrending transformation, with covariates
res <- lwdid(walmart, "log_wholesale_emp", "cid", "year", gvar = "first_year",
             rolling = "detrend", method = "ipwra",
             controls = c("x1", "x2", "x3"), reps = 999, seed = 1)
print(res)
plot(res)          # event-study plot with sup-t bands (requires ggplot2)
```

Large-N options: `method`, `pre` (pre-periods used in the transform),
`never` (never-treated controls only), `attgt` (return ATT(g,t) cells),
`ydot` (return transformed outcomes), `reps`/`seed` (bootstrap),
`level`, `nose`.

**Scope.** The small-N path (default; Lee & Wooldridge 2026b) and the
large-N path (`method=`; Lee & Wooldridge 2026a) are both implemented.
The large-N estimators (RA/IPW/IPWRA), WATT(r) aggregation, and
inference are validated against the Stata `lwdid` v2.4.2 output on
`lw_smoking` and `lw_walmart` (see `data-raw/stata-reference/`).

### Common-timing designs: `post` vs `dvar`

For common-timing designs (`gvar = NULL`), supply a binary `post`
indicator. `post` can be used two ways, and
[`lwdid()`](https://xiangao.github.io/lwdidR/reference/lwdid.md)
distinguishes them:

- If `post` is the treatment-on indicator D_it (1 only for treated units
  in post periods), the ever-post units are the treated group — no extra
  argument needed (backwards-compatible behavior).
- If `post` is a *calendar* post indicator (1 for **all** units once t ≥
  first treated period), the treatment group lives in a separate column.
  Pass it via the new `dvar` argument,
  e.g. `lwdid(..., post = "post", dvar = "treated")`. If `dvar` is
  omitted,
  [`lwdid()`](https://xiangao.github.io/lwdidR/reference/lwdid.md)
  auto-detects an unambiguous unit-invariant `treat`/`treated`/`D`/`d`
  column and errors if it cannot resolve one.

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

## Replication: Lee & Wooldridge (2026b), Section 7.2

Castle Doctrine laws and log homicide rates:

| Method              | ATT    | SE     | t-stat |
|---------------------|--------|--------|--------|
| Demeaning (OLS SE)  | 0.0917 | 0.0571 | 1.607  |
| Demeaning (HC3 SE)  | 0.0917 | 0.0612 | 1.500  |
| Detrending (HC3 SE) | 0.0666 | 0.0550 | 1.210  |

All results match the Castle Doctrine estimates reported in the text of
Section 7.2 of Lee & Wooldridge (2026b) within tolerance 0.001 (data set
used by Cunningham 2021). See
[`vignette("castle_law")`](https://xiangao.github.io/lwdidR/articles/)
for full replication.

## Main functions

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

## Citation

If you use `lwdidR`, please cite the Lee & Wooldridge papers and
acknowledge the original Stata package.

> Lee, S. J., & Wooldridge, J. M. (2026a). Simple Transformation
> Approach to Difference-in-Differences Estimation for Panel Data.
> *Journal of Business & Economic Statistics*, 1–27.
> <https://doi.org/10.1080/07350015.2026.2683047>
>
> Lee, S. J., & Wooldridge, J. M. (2026b). Simple Approaches to
> Inference with Difference-in-Differences Estimators with Small
> Cross-Sectional Sample Sizes. Working paper, SSRN 5325686.
> <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686>

**Original Stata package:** `lwdid` by Soo Jeong Lee and Jeffrey M.
Wooldridge, <https://github.com/Soo-econ/lwdid>. `lwdidR` is an
independent R port of that command (v1.0, based on 2026b). Section and
equation numbers cited for the Castle Doctrine example refer to Lee &
Wooldridge (2026b).
