# Castle Law Replication: Lee & Wooldridge (2025) Section 7.2

## Overview

This vignette replicates the staggered DiD estimates from Lee &
Wooldridge (2025), Section 7.2, using the `lwdidR` package. The analysis
studies the effect of Castle Doctrine laws on log homicide rates for US
states, 2000–2010.

**Expected results (from paper and Python replication):**

| Estimator           | ATT     | SE      | t-stat |
|---------------------|---------|---------|--------|
| Demeaning (OLS SE)  | ≈ 0.092 | ≈ 0.057 | ≈ 1.61 |
| Demeaning (HC3 SE)  | ≈ 0.092 | ≈ 0.061 | ≈ 1.50 |
| Detrending (HC3 SE) | ≈ 0.067 | ≈ 0.055 | ≈ 1.21 |

------------------------------------------------------------------------

## Setup

``` r

library(lwdidR)
library(knitr)
```

``` r

data_path <- system.file("extdata", "castle.csv", package = "lwdidR")
stopifnot(nchar(data_path) > 0)
castle <- read.csv(data_path)

# Construct gvar: first treatment year (NA for never treated)
castle$gvar <- castle$effyear
castle$gvar[is.na(castle$gvar) | castle$gvar == 0] <- NA

# Verify data structure
cat("Dimensions:", nrow(castle), "x", ncol(castle), "\n")
```

    ## Dimensions: 550 x 61

``` r

cat("Years:", min(castle$year), "to", max(castle$year), "\n")
```

    ## Years: 2000 to 2010

``` r

cat("States:", length(unique(castle$sid)), "\n")
```

    ## States: 50

``` r

cat("Treatment cohorts:\n")
```

    ## Treatment cohorts:

``` r

cohort_counts <- table(castle$gvar[!duplicated(castle$sid)], useNA = "ifany")
print(cohort_counts)
```

    ## 
    ## 2005 2006 2007 2008 2009 <NA> 
    ##    1   13    4    2    1   29

------------------------------------------------------------------------

## Demeaning Estimator

### With OLS Standard Errors

``` r

res_demean_ols <- lwdid(
  data          = castle,
  y             = "lhomicide",
  ivar          = "sid",
  tvar          = "year",
  gvar          = "gvar",
  rolling       = "demean",
  control_group = "never_treated",
  aggregate     = "overall",
  vce           = NULL  # homoskedastic OLS
)
print(res_demean_ols)
```

    ## 
    ## Lee-Wooldridge DiD (lwdidR)
    ## Design:      staggered
    ## Transf.:     demean
    ## VCE:         OLS
    ## --------------------------------------------------
    ## Overall ATT:   0.0917
    ## SE:            0.0571
    ## t-stat:        1.6067
    ## p-value:       0.1147
    ## --------------------------------------------------
    ## 
    ## Cohort-specific effects:
    ##  cohort     att      se  tstat pvalue
    ##    2005 0.08017 0.17305 0.4632 0.6468
    ##    2006 0.06824 0.07220 0.9450 0.3503
    ##    2007 0.11406 0.08998 1.2676 0.2144
    ##    2008 0.14605 0.13963 1.0459 0.3042
    ##    2009 0.21108 0.19105 1.1049 0.2786

### With HC3 Standard Errors

``` r

res_demean_hc3 <- lwdid(
  data          = castle,
  y             = "lhomicide",
  ivar          = "sid",
  tvar          = "year",
  gvar          = "gvar",
  rolling       = "demean",
  control_group = "never_treated",
  aggregate     = "overall",
  vce           = "hc3"
)
print(res_demean_hc3)
```

    ## 
    ## Lee-Wooldridge DiD (lwdidR)
    ## Design:      staggered
    ## Transf.:     demean
    ## VCE:         HC3
    ## --------------------------------------------------
    ## Overall ATT:   0.0917
    ## SE:            0.0612
    ## t-stat:        1.4997
    ## p-value:       0.1402
    ## --------------------------------------------------
    ## 
    ## Cohort-specific effects:
    ##  cohort     att      se tstat    pvalue
    ##    2005 0.08017 0.03215 2.493 1.884e-02
    ##    2006 0.06824 0.08920 0.765 4.488e-01
    ##    2007 0.11406 0.09838 1.159 2.552e-01
    ##    2008 0.14605 0.08203 1.780 8.548e-02
    ##    2009 0.21108 0.03550 5.946 2.115e-06

------------------------------------------------------------------------

## Detrending Estimator

``` r

res_detrend_hc3 <- lwdid(
  data          = castle,
  y             = "lhomicide",
  ivar          = "sid",
  tvar          = "year",
  gvar          = "gvar",
  rolling       = "detrend",
  control_group = "never_treated",
  aggregate     = "overall",
  vce           = "hc3"
)
print(res_detrend_hc3)
```

    ## 
    ## Lee-Wooldridge DiD (lwdidR)
    ## Design:      staggered
    ## Transf.:     detrend
    ## VCE:         HC3
    ## --------------------------------------------------
    ## Overall ATT:   0.0666
    ## SE:            0.0550
    ## t-stat:        1.2102
    ## p-value:       0.2321
    ## --------------------------------------------------
    ## 
    ## Cohort-specific effects:
    ##  cohort       att      se    tstat   pvalue
    ##    2005  0.139526 0.06496  2.14797 0.040513
    ##    2006  0.107340 0.05758  1.86411 0.069657
    ##    2007 -0.002499 0.14025 -0.01782 0.985897
    ##    2008 -0.126735 0.13892 -0.91231 0.369129
    ##    2009  0.126083 0.04250  2.96646 0.006102

------------------------------------------------------------------------

## Summary Comparison Table

``` r

paper_vals <- data.frame(
  Method       = c("Demeaning (OLS)", "Demeaning (HC3)", "Detrending (HC3)"),
  Paper_ATT    = c(0.092, 0.092, 0.067),
  Paper_SE     = c(0.057, NA, 0.055),
  Paper_tstat  = c(1.61, 1.50, 1.21)
)

r_vals <- data.frame(
  Method     = c("Demeaning (OLS)", "Demeaning (HC3)", "Detrending (HC3)"),
  R_ATT      = round(c(res_demean_ols$att_overall,
                        res_demean_hc3$att_overall,
                        res_detrend_hc3$att_overall), 4),
  R_SE       = round(c(res_demean_ols$se_overall,
                        res_demean_hc3$se_overall,
                        res_detrend_hc3$se_overall), 4),
  R_tstat    = round(c(res_demean_ols$tstat,
                        res_demean_hc3$tstat,
                        res_detrend_hc3$tstat), 3)
)

results_table <- merge(paper_vals, r_vals, by = "Method")
kable(results_table, caption = "Table 7.2 Replication: Lee & Wooldridge (2025)")
```

| Method           | Paper_ATT | Paper_SE | Paper_tstat |  R_ATT |   R_SE | R_tstat |
|:-----------------|----------:|---------:|------------:|-------:|-------:|--------:|
| Demeaning (HC3)  |     0.092 |       NA |        1.50 | 0.0917 | 0.0612 |   1.500 |
| Demeaning (OLS)  |     0.092 |    0.057 |        1.61 | 0.0917 | 0.0571 |   1.607 |
| Detrending (HC3) |     0.067 |    0.055 |        1.21 | 0.0666 | 0.0550 |   1.210 |

Table 7.2 Replication: Lee & Wooldridge (2025) {.table}

------------------------------------------------------------------------

## Cohort-Specific Effects

``` r

res_cohort <- lwdid(
  data          = castle,
  y             = "lhomicide",
  ivar          = "sid",
  tvar          = "year",
  gvar          = "gvar",
  rolling       = "demean",
  control_group = "never_treated",
  aggregate     = "cohort",
  vce           = "hc3"
)

cat("\nCohort-specific ATT estimates:\n")
```

    ## 
    ## Cohort-specific ATT estimates:

``` r

kable(
  res_cohort$att_by_cohort[, c("cohort", "att", "se", "ci_lower", "ci_upper",
                                 "tstat", "pvalue", "n_periods", "n_units")],
  digits = 4,
  caption = "Cohort Effects (Equation 7.9, Lee & Wooldridge 2025)"
)
```

|      | cohort |    att |     se | ci_lower | ci_upper |  tstat | pvalue | n_periods | n_units |
|:-----|-------:|-------:|-------:|---------:|---------:|-------:|-------:|----------:|--------:|
| 2005 |   2005 | 0.0802 | 0.0322 |   0.0143 |   0.1460 | 2.4932 | 0.0188 |         6 |       1 |
| 2006 |   2006 | 0.0682 | 0.0892 |  -0.1120 |   0.2485 | 0.7650 | 0.4488 |         5 |      13 |
| 2007 |   2007 | 0.1141 | 0.0984 |  -0.0866 |   0.3147 | 1.1594 | 0.2552 |         4 |       4 |
| 2008 |   2008 | 0.1460 | 0.0820 |  -0.0217 |   0.3138 | 1.7805 | 0.0855 |         3 |       2 |
| 2009 |   2009 | 0.2111 | 0.0355 |   0.1384 |   0.2838 | 5.9463 | 0.0000 |         2 |       1 |

Cohort Effects (Equation 7.9, Lee & Wooldridge 2025) {.table
style="width:100%;"}

------------------------------------------------------------------------

## All (g, r) Effects

``` r

res_none <- lwdid(
  data          = castle,
  y             = "lhomicide",
  ivar          = "sid",
  tvar          = "year",
  gvar          = "gvar",
  rolling       = "demean",
  control_group = "never_treated",
  aggregate     = "none",
  vce           = "hc3"
)

cat(sprintf("\nTotal (g,r) pairs estimated: %d\n", nrow(res_none$att_by_cohort_time)))
```

    ## 
    ## Total (g,r) pairs estimated: 20

``` r

kable(
  head(res_none$att_by_cohort_time[, c("cohort", "period", "event_time",
                                         "att", "se", "pvalue")], 15),
  digits = 4,
  caption = "(g, r)-Specific Effects (Equation 7.8, Lee & Wooldridge 2025)"
)
```

| cohort | period | event_time |     att |     se | pvalue |
|-------:|-------:|-----------:|--------:|-------:|-------:|
|   2005 |   2005 |          0 | -0.1332 | 0.0283 | 0.0001 |
|   2005 |   2006 |          1 |  0.0861 | 0.0345 | 0.0187 |
|   2005 |   2007 |          2 |  0.1640 | 0.0469 | 0.0016 |
|   2005 |   2008 |          3 |  0.1367 | 0.0507 | 0.0117 |
|   2005 |   2009 |          4 |  0.1284 | 0.0448 | 0.0079 |
|   2005 |   2010 |          5 |  0.0990 | 0.0488 | 0.0520 |
|   2006 |   2006 |          0 |  0.0663 | 0.0839 | 0.4342 |
|   2006 |   2007 |          1 |  0.1186 | 0.0935 | 0.2119 |
|   2006 |   2008 |          2 |  0.0220 | 0.1141 | 0.8478 |
|   2006 |   2009 |          3 |  0.0871 | 0.1077 | 0.4230 |
|   2006 |   2010 |          4 |  0.0471 | 0.0838 | 0.5771 |
|   2007 |   2007 |          0 |  0.1311 | 0.1971 | 0.5110 |
|   2007 |   2008 |          1 | -0.0767 | 0.1511 | 0.6152 |
|   2007 |   2009 |          2 |  0.2567 | 0.1506 | 0.0983 |
|   2007 |   2010 |          3 |  0.1452 | 0.1460 | 0.3277 |

(g, r)-Specific Effects (Equation 7.8, Lee & Wooldridge 2025) {.table}

------------------------------------------------------------------------

## Numerical Checks

``` r

cat("=== Numerical Verification ===\n\n")
```

    ## === Numerical Verification ===

``` r

tol <- 0.001
checks <- list(
  list(label = "Demeaning ATT",   got = res_demean_ols$att_overall, ref = 0.092),
  list(label = "OLS SE",          got = res_demean_ols$se_overall,  ref = 0.057),
  list(label = "OLS t-stat",      got = res_demean_ols$tstat,       ref = 1.607),
  list(label = "HC3 t-stat",      got = res_demean_hc3$tstat,       ref = 1.50),
  list(label = "Detrending ATT",  got = res_detrend_hc3$att_overall,ref = 0.067),
  list(label = "Detrend HC3 t",   got = res_detrend_hc3$tstat,      ref = 1.21)
)

all_pass <- TRUE
for (chk in checks) {
  diff  <- abs(chk$got - chk$ref)
  pass  <- diff < tol
  all_pass <- all_pass && pass
  cat(sprintf("%-20s  got: %7.4f  ref: %5.3f  diff: %.5f  %s\n",
              chk$label, chk$got, chk$ref, diff,
              if (pass) "PASS" else "FAIL"))
}
```

    ## Demeaning ATT         got:  0.0917  ref: 0.092  diff: 0.00025  PASS
    ## OLS SE                got:  0.0571  ref: 0.057  diff: 0.00010  PASS
    ## OLS t-stat            got:  1.6067  ref: 1.607  diff: 0.00033  PASS
    ## HC3 t-stat            got:  1.4997  ref: 1.500  diff: 0.00026  PASS
    ## Detrending ATT        got:  0.0666  ref: 0.067  diff: 0.00045  PASS
    ## Detrend HC3 t         got:  1.2102  ref: 1.210  diff: 0.00024  PASS

``` r

cat(sprintf("\nAll checks pass: %s\n", all_pass))
```

    ## 
    ## All checks pass: TRUE

------------------------------------------------------------------------

## Session Info

``` r

sessionInfo()
```

    ## R version 4.6.0 (2026-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.51   lwdidR_0.1.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.57         cachem_1.1.0      htmltools_0.5.9   rmarkdown_2.31   
    ##  [9] lifecycle_1.0.5   cli_3.6.6         sass_0.4.10       pkgdown_2.2.0    
    ## [13] textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2 compiler_4.6.0   
    ## [17] tools_4.6.0       ragg_1.5.2        evaluate_1.0.5    bslib_0.11.0     
    ## [21] yaml_2.3.12       jsonlite_2.0.0    rlang_1.2.0       fs_2.1.0
