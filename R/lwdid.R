# Main dispatch function for lwdidR

#' Lee-Wooldridge DiD via unit-specific pre-treatment transformations
#'
#' Main entry point implementing the Lee & Wooldridge (2026) panel
#' difference-in-differences estimator. Supports both common-timing and
#' staggered adoption designs. Each unit's outcome is residualised using
#' only its own pre-treatment observations before running a pooled
#' cross-sectional OLS.
#'
#' @param data A long-format panel data frame (one row per unit-period).
#' @param y Character. Name of the outcome column.
#' @param ivar Character. Name of the unit identifier column.
#' @param tvar Character. Name of the calendar time column (numeric or integer).
#' @param gvar Character or NULL. Name of the first-treatment-year column for
#'   staggered designs. Units with value `0` or `NA` are treated as never-treated.
#'   Set to `NULL` (default) for common-timing designs (supply `post` instead).
#' @param post Character or NULL. Name of a binary post-treatment indicator column
#'   (0 = pre, 1 = post). Required when `gvar = NULL`; ignored otherwise.
#' @param dvar Character or NULL. Name of a unit-level treatment-group indicator
#'   for common-timing designs. If omitted, `lwdid()` keeps backwards-compatible
#'   behavior when `post` is the treatment-on indicator; if `post` is a calendar
#'   post indicator, it looks for an unambiguous unit-invariant treatment column
#'   named `treat`, `treated`, `D`, or `d`.
#' @param rolling Character. Transformation method applied to each unit's
#'   pre-treatment observations:
#'   \describe{
#'     \item{`"demean"`}{Subtract the unit's pre-period mean (default).}
#'     \item{`"detrend"`}{Remove a linear trend fitted on pre-periods.}
#'     \item{`"demeanq"`}{Seasonal demeaning; requires `season_var`.}
#'     \item{`"detrendq"`}{Seasonal detrending; requires `season_var`.}
#'   }
#' @param control_group Character. Control group for staggered designs:
#'   `"never_treated"` (default) or `"not_yet_treated"`.
#' @param aggregate Character. Aggregation level for staggered designs:
#'   `"overall"` (default), `"cohort"`, or `"none"` (returns all (g,r) pairs).
#' @param vce Character or NULL. Variance-covariance estimator:
#'   `NULL` (homoskedastic OLS), `"hc1"`, `"hc3"`, `"cluster"`,
#'   `"wildboot"` (wild cluster bootstrap), or `"permutation"` (randomisation
#'   inference). The last two are distribution-free and recommended at small N.
#' @param cluster_var Character or NULL. Column name for clustering; required
#'   when `vce = "cluster"`.
#' @param nboot Integer. Number of bootstrap replications for `vce = "wildboot"`
#'   (default 999).
#' @param nperm Integer. Number of permutations for `vce = "permutation"`
#'   (default 999).
#' @param vce_inner Character. Inner variance estimator used when computing the
#'   observed t-statistic inside the wild bootstrap (default `"hc3"`).
#' @param controls Character vector or NULL. Names of time-invariant control
#'   variables to include in the cross-sectional regression.
#' @param season_var Character or NULL. Column name of the seasonal indicator
#'   (required for `rolling = "demeanq"` or `"detrendq"`).
#' @param method Character or NULL. If one of `"ra"`, `"ipw"`, `"ipwra"`, the
#'   large-N path (Lee & Wooldridge 2026a) is used: per-cohort/period ATT(g,t)
#'   via regression adjustment, inverse-probability weighting, or doubly-robust
#'   IPWRA, aggregated to the event-study WATT(r) path with wild-cluster-bootstrap
#'   inference. `ipw`/`ipwra` require `controls`. `rolling` must be `"demean"` or
#'   `"detrend"`. If `NULL` (default), the small-N path is used.
#' @param pre Integer. Large-N only. Number of most-recent pre-periods used in the
#'   transformation; `-1` (default) uses all pre-periods.
#' @param never Logical. Large-N only. If `TRUE`, the comparison group is
#'   never-treated units only (default `FALSE` uses not-yet-treated units too).
#' @param attgt Logical. Large-N only. If `TRUE`, the ATT(g,t) cell estimates are
#'   returned (and printed).
#' @param ydot Logical. Large-N only. If `TRUE`, the per-cohort transformed
#'   outcomes are returned.
#' @param reps Integer. Large-N only. Wild cluster bootstrap replications (default 999).
#' @param seed Integer or NULL. Large-N only. Seed for the wild bootstrap.
#' @param level Numeric. Confidence level for the large-N path (default 95).
#' @param nose Logical. Large-N only. If `TRUE`, skip standard errors (faster).
#'
#' @return An object of class `"lwdid"`, a list containing:
#'   \describe{
#'     \item{`design`}{`"staggered"` or `"common_timing"`.}
#'     \item{`att_overall`}{Estimated overall ATT.}
#'     \item{`se_overall`}{Standard error of overall ATT.}
#'     \item{`tstat`}{t-statistic.}
#'     \item{`pvalue`}{Two-sided p-value.}
#'     \item{`att_by_cohort`}{Data frame of cohort-specific ATTs (staggered only).}
#'     \item{`att_by_cohort_time`}{Data frame of (g,r)-specific ATTs (staggered only).}
#'     \item{`att_by_period`}{Data frame of period-specific ATTs (common timing only).}
#'     \item{`ci_lower`, `ci_upper`}{95% confidence interval bounds (common timing only).}
#'     \item{`N`}{Sample size at first post-treatment period (common timing only).}
#'   }
#'
#' @references
#'   Lee, S. J., & Wooldridge, J. M. (2026a). Simple Transformation Approach to
#'   Difference-in-Differences Estimation for Panel Data. Journal of Business &
#'   Economic Statistics, 1-27. \doi{10.1080/07350015.2026.2683047}
#'
#'   Lee, S. J., & Wooldridge, J. M. (2026b). Simple Approaches to Inference with
#'   Difference-in-Differences Estimators with Small Cross-Sectional Sample Sizes.
#'   Working paper, SSRN 5325686.
#'   \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325686}
#'
#'   \code{lwdidR} is an R port of the original Stata \code{lwdid} command by
#'   Soo Jeong Lee and Jeffrey M. Wooldridge:
#'   \url{https://github.com/Soo-econ/lwdid}. The Castle Doctrine example (Section
#'   7.2) and the small-sample inference options are from (2026b).
#'
#' @examples
#' # Load bundled Castle Doctrine dataset
#' castle <- read.csv(system.file("extdata", "castle.csv", package = "lwdidR"))
#' castle$gvar <- castle$effyear
#' castle$gvar[is.na(castle$gvar) | castle$gvar == 0] <- NA
#'
#' # Staggered design with demeaning and HC3 standard errors
#' res <- lwdid(castle, "lhomicide", "sid", "year",
#'              gvar = "gvar", rolling = "demean", vce = "hc3")
#' print(res)
#'
#' @export
lwdid <- function(data, y, ivar, tvar,
                  gvar = NULL, post = NULL, dvar = NULL,
                  rolling = "demean",
                  method = NULL,
                  control_group = "never_treated",
                  aggregate = "overall",
                  vce = NULL, cluster_var = NULL,
                  controls = NULL,
                  season_var = NULL,
                  pre = -1L, never = FALSE, attgt = FALSE, ydot = FALSE,
                  reps = 999L, seed = NULL, level = 95, nose = FALSE,
                  nboot = 999, nperm = 999, vce_inner = "hc3") {

  # Input validation
  stopifnot(is.data.frame(data))
  for (v in c(y, ivar, tvar)) {
    if (!v %in% names(data)) stop(sprintf("Column '%s' not found in data.", v))
  }

  # --- Large-N path (Lee & Wooldridge 2026a): triggered by method = ra/ipw/ipwra
  if (!is.null(method)) {
    if (is.null(gvar)) stop("Large-N path (method=) requires 'gvar'.")
    if (!gvar %in% names(data)) stop(sprintf("Column '%s' not found in data.", gvar))
    if (!method %in% c("ra", "ipw", "ipwra"))
      stop("method must be 'ra', 'ipw', or 'ipwra'.")
    return(.lwdid_large_dispatch(
      data = data, y = y, ivar = ivar, tvar = tvar, gvar = gvar,
      rolling = rolling, method = method, controls = controls,
      pre = pre, never = never, attgt = attgt, ydot = ydot,
      cluster_var = cluster_var, reps = reps, seed = seed,
      level = level, nose = nose))
  }
  if (!rolling %in% c("demean", "detrend", "demeanq", "detrendq")) {
    stop("rolling must be one of: demean, detrend, demeanq, detrendq")
  }
  valid_vce_all <- c("hc1", "hc3", "cluster", "wildboot", "permutation")
  if (!is.null(vce) && !tolower(vce) %in% valid_vce_all) {
    stop(sprintf("vce must be NULL or one of: %s", paste(valid_vce_all, collapse = ", ")))
  }
  if (!is.null(vce) && tolower(vce) == "cluster") {
    if (is.null(cluster_var)) stop("vce='cluster' requires cluster_var.")
    if (!cluster_var %in% names(data)) stop(sprintf("Column '%s' not found in data.", cluster_var))
  }

  # Dispatch
  if (!is.null(gvar)) {
    # --- Staggered adoption path ---
    if (!gvar %in% names(data)) stop(sprintf("Column '%s' not found in data.", gvar))

    res <- lwdid_staggered(
      data          = data,
      y             = y,
      ivar          = ivar,
      tvar          = tvar,
      gvar          = gvar,
      rolling       = rolling,
      control_group = control_group,
      aggregate     = aggregate,
      vce           = vce,
      cluster_var   = cluster_var,
      controls      = controls,
      season_var    = season_var,
      nboot         = nboot,
      nperm         = nperm,
      vce_inner     = vce_inner
    )

    out <- structure(
      list(
        design          = "staggered",
        att_overall     = res$att_overall$att,
        se_overall      = res$att_overall$se,
        tstat           = res$att_overall$tstat,
        pvalue          = res$att_overall$pvalue,
        att_by_cohort   = res$att_by_cohort,
        att_by_cohort_time = res$att_by_cohort_time,
        cohorts         = res$cohorts,
        rolling         = rolling,
        control_group   = control_group,
        aggregate       = aggregate,
        vce             = vce,
        y = y, ivar = ivar, tvar = tvar, gvar = gvar
      ),
      class = "lwdid"
    )

  } else {
    # --- Common timing path ---
    if (is.null(post)) stop("For common-timing design, supply 'post' column name.")
    if (!post %in% names(data)) stop(sprintf("Column '%s' not found in data.", post))

    all_periods  <- sort(unique(data[[tvar]]))
    post_periods <- sort(unique(data[[tvar]][data[[post]] == 1]))
    if (length(post_periods) == 0) stop("No post-treatment periods found.")
    tpost1       <- min(post_periods)

    # Treatment-group indicator. Older callers sometimes used `post` as the
    # treatment-on indicator D_it; newer/common DiD data usually use `post` as a
    # calendar post indicator and keep treatment status in a separate column.
    if (!is.null(dvar)) {
      if (!dvar %in% names(data)) stop(sprintf("Column '%s' not found in data.", dvar))
      treated_units <- unique(data[[ivar]][data[[dvar]] == 1])
    } else {
      post_by_time <- stats::aggregate(data[[post]], list(data[[tvar]]), function(x) {
        ux <- unique(stats::na.omit(x))
        if (length(ux) == 1) ux else NA
      })[[2]]
      post_is_calendar <- all(!is.na(post_by_time)) &&
        any(post_by_time == 0) && any(post_by_time == 1)

      if (post_is_calendar) {
        candidate_names <- intersect(c("treat", "treated", "D", "d"), names(data))
        candidate_names <- candidate_names[vapply(candidate_names, function(v) {
          vals <- unique(stats::na.omit(data[[v]]))
          unit_vals <- stats::aggregate(data[[v]], list(data[[ivar]]), function(x) {
            ux <- unique(stats::na.omit(x))
            length(ux) == 1
          })[[2]]
          length(vals) == 2 && all(sort(vals) == c(0, 1)) && all(unit_vals)
        }, logical(1))]
        if (length(candidate_names) != 1) {
          stop("For common-timing designs with a calendar post indicator, supply 'dvar' for the treatment-group indicator.")
        }
        dvar <- candidate_names[[1]]
        treated_units <- unique(data[[ivar]][data[[dvar]] == 1])
      } else {
        # Backwards-compatible path: `post` is D_it, so ever-post units are treated.
        treated_units <- unique(data[[ivar]][data[[post]] == 1])
      }
    }

    # Calendar post indicator: 1 for ALL units once t >= tpost1.
    # This is passed to apply_transform so that control units (which always
    # have post == 0 in the user-supplied column) still get ydot_postavg
    # computed and appear in the firstpost cross-section.
    data$.post_cal_ <- as.integer(data[[tvar]] >= tpost1)

    df_trans <- apply_transform(
      df        = data,
      y         = y,
      ivar      = ivar,
      tindex    = tvar,
      post      = ".post_cal_",
      rolling   = rolling,
      tpost1    = tpost1,
      season_var = season_var
    )

    df_trans$d_ <- as.integer(df_trans[[ivar]] %in% treated_units)

    att_res <- estimate_att(df_trans, d = "d_", vce = vce,
                            cluster_var = cluster_var, controls = controls,
                            nboot = nboot, nperm = nperm, vce_inner = vce_inner)

    per_res <- estimate_period_effects(
      df          = df_trans,
      d           = "d_",
      tindex      = tvar,
      post_periods = post_periods,
      vce         = vce,
      cluster_var  = cluster_var,
      controls    = controls,
      nboot       = nboot,
      nperm       = nperm,
      vce_inner   = vce_inner
    )

    out <- structure(
      list(
        design       = "common_timing",
        att_overall  = att_res$att,
        se_overall   = att_res$se,
        tstat        = att_res$tstat,
        pvalue       = att_res$pvalue,
        ci_lower     = att_res$ci_lower,
        ci_upper     = att_res$ci_upper,
        att_by_period = per_res,
        N            = att_res$N,
        rolling      = rolling,
        vce          = vce,
        y = y, ivar = ivar, tvar = tvar, post = post, dvar = dvar
      ),
      class = "lwdid"
    )
  }

  out
}
