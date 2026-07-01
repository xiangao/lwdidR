# lwdidR — project notes for Claude

R implementation of the Lee & Wooldridge panel DiD estimator: residualize each
unit's outcome using only its *own* pre-treatment observations (demean /
detrend / demeanq / detrendq), then run pooled cross-sectional OLS on the
residualized outcome. It is an R port of the original Stata `lwdid` command
(Soo Jeong Lee & Jeffrey M. Wooldridge, <https://github.com/Soo-econ/lwdid>).

Two estimator families live behind one entry point, `lwdid()`:

- **Small-N path** (default, `method = NULL`) — Lee & Wooldridge (2026b,
  SSRN 5325686). Pooled cross-sectional OLS on `ydot_postavg`, with
  homoskedastic/HC1/HC3/cluster SEs or distribution-free wild-cluster-bootstrap
  / permutation inference. This is the whole original package (v0.1.x).
  Implemented in `R/lwdid.R`, `R/transform.R`, `R/estimate.R`, `R/staggered.R`,
  `R/methods.R`.
- **Large-N path** (`method = "ra"`/`"ipw"`/`"ipwra"`) — Lee & Wooldridge
  (2026a, JBES, doi:10.1080/07350015.2026.2683047). Added in v0.2.0. Estimates
  ATT(g,t) per cohort × calendar period via regression adjustment / IPW /
  doubly-robust IPWRA, aggregates (treated-count-weighted) to an event-study
  WATT(r) path plus Pre/Post averages, with wild-cluster-bootstrap SEs and
  sup-t simultaneous confidence bands. Implemented in `R/large_n.R`,
  `R/transform_large.R`, `R/estimators_large.R`, `R/aggregate_large.R`,
  `R/inference_large.R`, `R/methods_large.R`, `R/plot_large.R`. Result class
  is `c("lwdid_large", "lwdid")`.

Both paths are ported line-by-line from Stata `lwdid` v2.4.2's `.ado`/Mata
source (see file-header comments in each `R/*_large.R` file for the
approximate source line ranges) — this is not a reimplementation from the
papers alone, it is a faithful transcription, so when in doubt about intended
behavior, the fidelity target is Stata's code, not a plausible-looking
rewrite.

## Running tests

```r
devtools::load_all(quiet = TRUE)
testthat::test_dir("tests/testthat")
```

Current status: **65/65 pass, 0 skipped** (confirmed 2026-07-01; runtime is a
minute or two because of bootstrap-based inference tests: wild-cluster
bootstrap with `reps ~ 99-999` runs inside several tests).

The large-N Stata cross-validation tests (`test-large_n.R`) read
`data-raw/stata-reference/ref_walmart_*.csv` and `ref_smoking_smallN.csv` if
present, and `skip()` cleanly if absent (e.g. on CI without Stata). They were
present and passing in this run — do not assume they always run; check for
`SKIP` in the summary line if regenerating them isn't possible.

## Stata cross-validation (`data-raw/stata-reference/`)

- `lwdid_reference.do` regenerates the oracle: installs Stata `lwdid`, fetches
  its bundled `lw_smoking.dta` and `lw_walmart.dta`, and writes 9 large-N
  configs (`ref_walmart_<config>.csv`, covering demean/detrend ×
  ra/ipw/ipwra, plus `_nox` no-controls, `_pre3` limited pre-window, `_never`
  never-treated-only variants) and one small-N config
  (`ref_smoking_smallN.csv`).
- All generated CSVs are git-ignored (regenerate via the do-file if Stata is
  available; `lw_walmart` also ships pre-baked at `inst/extdata/lw_walmart.rds`
  since it belongs to the upstream `lwdid` package, not this one).
- Validated tolerance (per `data-raw/stata-reference/README.md`, Stata
  `lwdid` v2.4 vs `lwdidR` 0.2.0): WATT(r)/Pre/Post point estimates match to
  <= 1.5e-8 (deterministic — no RNG involved in point estimates); bootstrap SEs
  agree to within ~3% median relative gap (not exact, since Stata's and R's
  bootstrap RNGs differ — see below).
- `test-lwdid.R` has no Stata oracle; it validates against known DGP truths
  (synthetic panels with a specified ATT) and structural invariants (SE > 0,
  correct S3 class, cohort tables non-empty, etc.), plus the Castle Doctrine
  replication numbers reported in the README/vignette.

### Why bootstrap SEs can't match Stata bit-for-bit (and what does)

Both wild-cluster-bootstrap implementations (small-N `.wild_bootstrap()` in
`R/estimate.R`, large-N `.inference_large()` in `R/inference_large.R`) draw
independent Rademacher weights, so no seed alignment will make R and Stata
agree exactly. Instead, `.inference_large()` also computes a **deterministic**
analytic cluster-robust SE (`se_analytic` = sqrt of the cluster-sum-of-squares
of the centered influence function) as the true large-sample limit of the
Rademacher bootstrap variance, and the Stata cross-validation test compares
that analytic SE to Stata's reported bootstrap SE (median relative gap <
6%) — this is the correct invariant to check, not bootstrap-SE-vs-bootstrap-SE
matching.

## The `post` vs `dvar` distinction (common-timing designs only; `gvar = NULL`)

This is the easiest thing to get backwards. `lwdid()` accepts a `post` column
for common-timing designs, and it can mean **two different things**:

1. **`post` = D_it, the treatment-ON indicator** (1 only for treated units in
   post periods, 0 everywhere else including all periods for control units).
   This is the original/backwards-compatible behavior: no `dvar` needed.
   `ever_post` units (any row with `post==1`) are inferred to be the treated
   group.
2. **`post` = a calendar indicator** (1 for ALL units, treated and control
   alike, once `t >= first treated period`; 0 before). In this case treatment
   status lives in a *separate* column, and you MUST either pass it via
   `dvar = "treated"` or have an unambiguous unit-invariant 0/1 column named
   exactly `treat`, `treated`, `D`, or `d` in the data (auto-detected).

`lwdid()` disambiguates automatically by checking, for each calendar period,
whether `post` takes a single value across ALL units — if some periods are
"mixed" (both 0s and 1s for the same `t`), `post` is being used as D_it
(case 1); if `post` is constant-by-period and takes both 0 and 1 across
different periods, it's a calendar indicator (case 2), and `dvar` resolution
kicks in (see `R/lwdid.R` lines ~204-242, `.post_cal_` internal column). If
you build a calendar-`post` dataset without an obvious `treat`-like column and
without passing `dvar`, `lwdid()` throws rather than guessing — this is
deliberate, do not "fix" it by picking a default.

Internally, regardless of which case applies, `lwdid()` builds its own
`.post_cal_` column (`t >= tpost1` for every unit) before calling
`apply_transform()`, because control units need `ydot_postavg` computed too
(their `post` column may be all-0 under case 1) so they can appear in the
`firstpost` cross-section used for `estimate_att()`.

## Other real gotchas found while reading `R/`

- **Anchoring is mechanical, not estimated.** In the large-N path, the
  transformed outcome is forced to exactly 0 at the normalization period(s):
  `r == -1` for demean, `r %in% c(-2, -1)` for detrend (`R/transform_large.R`).
  This propagates to `WATT(r)` (`.aggregate_large()` hard-codes `watt = 0` and
  `.inference_large()` hard-codes `se = 0` at those `r`) — a future change to
  one anchor list without the matching change in the other two files will
  silently break the sup-t band construction (the anchor's SE=0 column stays
  in `IF_r` unpermuted).
- **IPW-only excludes covariates from the outcome model.** For
  `method = "ipw"`, `.fit_cell_large()` sets `xc_out <- NULL` — covariates
  enter *only* through the logit propensity score, never as outcome-model
  regressors (mirrors Stata's `reg yvar dvar_g [aw=ipw]` with no `x` terms).
  For `ra`/`ipwra`, covariates are cohort-centered (`v - mean(v among
  first-treated-in-g)`, computed in `lwdid_large()`) before entering as `xc`
  and its interaction with `d`.
- **Control-group membership at cell (g,t) depends on whether t is pre or
  post the cohort's treatment**, and further on `never`
  (`R/large_n.R` lines ~66-72): pre-period cells (`tt < g`) include not-yet
  cohorts unconditionally (`gv > g`); post-period cells additionally require
  `gv > tt` (a cohort that already switched on by `tt` cannot serve as
  control). Getting this backwards silently contaminates the control group
  with already-treated units in post periods.
- **Overall ATT (small-N staggered) is not a naive average of cohort ATTs.**
  It is one pooled OLS: each treated unit contributes its own-cohort
  `ydot_postavg`; each never-treated unit contributes a treated-count-weighted
  average of `ydot_postavg` across *all* cohorts it could serve as control
  for. README explicitly flags this yields SE ≈ 0.061, not 0.051 from
  delta-method averaging — if a future refactor of `lwdid_staggered()`
  "simplifies" this into a two-stage average, it will quietly change the SE
  and break the Castle Doctrine replication numbers.
- **`prepare_controls()` silently drops controls, with only a `warning()`**,
  if `N1 <= K+1` or `N0 <= K+1` in the firstpost/period cross-section
  (`R/estimate.R`). A silently-underpowered cell in a period-by-period
  small-N regression can pass through without an error — check warnings when
  interpreting `att_by_period`/`att_by_cohort_time` on thin panels.
- **`detrend`/`detrendq` require >= 2 pre-periods per unit** (small-N path
  errors informatively; large-N path instead sets that unit's `ydot` to `NA`
  for the whole panel rather than erroring, so cohorts with a mix of short-
  and long-history units will silently lose the short-history units from
  estimation rather than failing).
- **`pre` (large-N only)** restricts the pre-window to the `pre` most recent
  pre-periods (`[g - pre, g - 1]`); `-1` (default) uses the full pre-history.
  This is a real Stata-parity option validated by the `_pre3` reference
  config, not dead code — don't assume "more pre-periods is always better"
  when tuning it; it trades off bias (a longer trend/level window) against
  variance and susceptibility to earlier confounding shocks.
- Field naming across the two paths is inconsistent by design (not a bug to
  "fix"): small-N staggered results use `att`/`se`/`ci_lower`/`ci_upper`;
  large-N results use `watt`/`se`/`lower_ci`/`upper_ci` plus `se_analytic`.
  They are genuinely different objects (`class = "lwdid"` vs
  `c("lwdid_large","lwdid")`) with different `print`/`summary` methods
  (`R/methods.R` vs `R/methods_large.R`) — don't try to unify the column
  names without checking every caller.

## Package layout

Not the workspace-standard `/code /data /output /report` layout — this is a
normal R package (`R/`, `tests/testthat/`, `man/`, `vignettes/`, `inst/`,
`data-raw/`). `data-raw/stata-reference/` holds the Stata oracle setup;
`inst/extdata/` ships the bundled `castle.csv` (Cunningham 2021, Section 7.2
replication) and `lw_walmart.rds`. Docs are pre-rendered and deployed via
pkgdown to <https://xiangao.github.io/lwdidR/> (see
`feedback_quarto_local_render_deploy`-style workflow: render vignettes
locally, commit, don't let CI re-render — check `.github/workflows/` before
assuming CI builds the site).

## Citation / attribution reminders

Two distinct papers, cite the right one for the right claim:

- Lee & Wooldridge (2026a), JBES — the outcome transformation itself
  (demean/detrend) and the large-N RA/IPW/IPWRA estimator.
- Lee & Wooldridge (2026b), SSRN 5325686 — small-sample inference (wild
  bootstrap, permutation) and the Castle Doctrine application (Section 7.2).

The DESCRIPTION/README got this wrong once already (v0.1.0 → v0.1.2 fixed a
misattributed author name and a wrong paper citation for the Castle Doctrine
numbers — see NEWS.md) — double check paper attribution before touching
docstrings or the README.
