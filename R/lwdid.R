# Main dispatch function for lwdidR

#' Lee-Wooldridge DiD via unit-specific pre-treatment transformations
#'
#' Main entry point implementing the Lee & Wooldridge (2025) panel
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
#'   `NULL` (homoskedastic OLS), `"hc1"`, `"hc3"`, or `"cluster"`.
#' @param cluster_var Character or NULL. Column name for clustering; required
#'   when `vce = "cluster"`.
#' @param controls Character vector or NULL. Names of time-invariant control
#'   variables to include in the cross-sectional regression.
#' @param season_var Character or NULL. Column name of the seasonal indicator
#'   (required for `rolling = "demeanq"` or `"detrendq"`).
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
#' @references Lee, Y., & Wooldridge, J. M. (2025). A simple panel data approach
#'   to difference-in-differences under general treatment patterns.
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
                  gvar = NULL, post = NULL,
                  rolling = "demean",
                  control_group = "never_treated",
                  aggregate = "overall",
                  vce = NULL, cluster_var = NULL,
                  controls = NULL,
                  season_var = NULL) {

  # Input validation
  stopifnot(is.data.frame(data))
  for (v in c(y, ivar, tvar)) {
    if (!v %in% names(data)) stop(sprintf("Column '%s' not found in data.", v))
  }
  if (!rolling %in% c("demean", "detrend", "demeanq", "detrendq")) {
    stop("rolling must be one of: demean, detrend, demeanq, detrendq")
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
      season_var    = season_var
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
    tpost1       <- min(post_periods)

    df_trans <- apply_transform(
      df        = data,
      y         = y,
      ivar      = ivar,
      tindex    = tvar,
      post      = post,
      rolling   = rolling,
      tpost1    = tpost1,
      season_var = season_var
    )

    # Treatment indicator: 1 if unit ever treated (post==1 for any period)
    treated_units <- unique(data[[ivar]][data[[post]] == 1])
    df_trans$d_ <- as.integer(df_trans[[ivar]] %in% treated_units)

    att_res <- estimate_att(df_trans, d = "d_", vce = vce,
                            cluster_var = cluster_var, controls = controls)

    per_res <- estimate_period_effects(
      df          = df_trans,
      d           = "d_",
      tindex      = tvar,
      post_periods = post_periods,
      vce         = vce,
      cluster_var  = cluster_var,
      controls    = controls
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
        y = y, ivar = ivar, tvar = tvar, post = post
      ),
      class = "lwdid"
    )
  }

  out
}
