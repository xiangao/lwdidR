# Staggered adoption LWDID estimator
# Mirrors Python lwdid/staggered/estimation.py + aggregation.py


#' LWDID for staggered adoption (gvar case)
#'
#' Python-equivalent algorithm:
#'   For each cohort g:
#'   1. Apply transformation using pre-periods (t < g) for ALL units
#'   2. Compute ydot_postavg(g,i) per unit
#'   Overall ATT: pooled regression
#'   - Treated unit i: Y_bar_i = ydot_postavg(g_i, i)
#'   - Never-treated unit i: Y_bar_i = sum_g w_g * ydot_postavg(g,i)
#'   - OLS: Y_bar ~ D (D=1 if ever treated) on N_treat + N_NT units
#'   Cohort ATT_g: per-cohort regression
#'   - Y_bar(g,i) = ydot_postavg(g,i) for cohort g units + NT controls
#'   - OLS: Y_bar_g ~ D_g
#'
#' @param data Long-format panel data
#' @param y Outcome column
#' @param ivar Unit identifier
#' @param tvar Calendar time column
#' @param gvar First treatment year (0/NA = never treated)
#' @param rolling Transformation method
#' @param control_group "never_treated" or "not_yet_treated"
#' @param aggregate "overall", "cohort", or "none"
#' @param vce Variance estimator
#' @param cluster_var Cluster variable
#' @param controls Control variable names
#' @param season_var Seasonal indicator column
#' @return Named list with att_overall, att_by_cohort, att_by_cohort_time
#' @keywords internal
lwdid_staggered <- function(data, y, ivar, tvar, gvar,
                             rolling       = "demean",
                             control_group = "never_treated",
                             aggregate     = "overall",
                             vce           = NULL,
                             cluster_var   = NULL,
                             controls      = NULL,
                             season_var    = NULL) {

  never_mask  <- is.na(data[[gvar]]) | data[[gvar]] == 0 | is.infinite(data[[gvar]])
  cohorts     <- sort(unique(data[[gvar]][!never_mask]))
  all_periods <- sort(unique(data[[tvar]]))
  Tmax        <- max(all_periods)

  # --- Step 1: compute ydot_postavg for each (cohort g, unit i) ---
  # ydot_g_postavg[unit_id] = mean(ydot | t >= g) using pre < g as pre-period
  # Stored as a named list indexed by cohort
  ydot_g_list  <- list()  # ydot_g_list[[as.char(g)]] = named vector, indexed by unit
  cohort_n1    <- integer(length(cohorts))  # N1 for each cohort
  names(cohort_n1) <- as.character(cohorts)

  for (g in cohorts) {
    pre_periods  <- all_periods[all_periods < g]
    post_periods <- all_periods[all_periods >= g]
    tpost1       <- min(post_periods)

    if (length(pre_periods) < 1) {
      warning(sprintf("Cohort %d: no pre-periods, skipping.", g)); next
    }
    if (rolling %in% c("detrend","detrendq") && length(pre_periods) < 2) {
      warning(sprintf("Cohort %d: need >=2 pre-periods for detrend, skipping.", g)); next
    }

    # Build cohort dataset
    treated_units <- unique(data[[ivar]][!never_mask & data[[gvar]] == g])

    if (control_group == "never_treated") {
      control_units <- unique(data[[ivar]][never_mask])
    } else {
      nyt_units     <- unique(data[[ivar]][!never_mask & data[[gvar]] > g])
      control_units <- unique(c(unique(data[[ivar]][never_mask]), nyt_units))
    }
    if (length(control_units) == 0) {
      warning(sprintf("Cohort %d: no controls, skipping.", g)); next
    }

    cohort_units <- c(treated_units, control_units)
    df_c <- data[data[[ivar]] %in% cohort_units, , drop = FALSE]
    df_c$post_g <- as.integer(df_c[[tvar]] >= g)
    df_c$d_g    <- as.integer(df_c[[ivar]] %in% treated_units)

    df_t <- apply_transform(df_c, y = y, ivar = ivar, tindex = tvar,
                            post = "post_g", rolling = rolling,
                            tpost1 = tpost1, season_var = season_var)

    # ydot_postavg per unit for this cohort
    # (apply_transform already computes ydot_postavg as mean of ydot | post_g=1)
    fp    <- df_t[df_t$firstpost, , drop = FALSE]
    ybar  <- stats::setNames(fp$ydot_postavg, fp[[ivar]])
    ydot_g_list[[as.character(g)]] <- ybar

    cohort_n1[as.character(g)] <- sum(fp$d_g == 1)
  }

  valid_cohorts <- cohorts[as.character(cohorts) %in% names(ydot_g_list)]
  if (length(valid_cohorts) == 0) stop("No cohorts could be estimated.")

  # --- Step 2: Cohort-level ATT_g ---
  # Per-cohort regression: ydot_postavg_g ~ D_g on (cohort g + NT controls)
  never_units   <- unique(data[[ivar]][never_mask])
  cohort_res    <- list()

  for (g in valid_cohorts) {
    gkey  <- as.character(g)
    ybar  <- ydot_g_list[[gkey]]
    units <- names(ybar)

    # Reconstruct firstpost sample
    treated_in_g <- units[units %in% unique(data[[ivar]][!never_mask & data[[gvar]] == g])]
    control_in_g <- units[units %in% never_units]
    sample_units <- c(treated_in_g, control_in_g)
    samp_ybar    <- ybar[sample_units]
    samp_d       <- as.integer(sample_units %in% treated_in_g)

    N1 <- sum(samp_d == 1); N0 <- sum(samp_d == 0)
    if (N1 == 0 || N0 == 0 || length(sample_units) < 3) next

    X   <- cbind(intercept = 1, d = samp_d)
    fit <- .ols_fit(X, samp_ybar)
    res <- .extract_inference(fit, vce,
              if (!is.null(vce) && tolower(vce) == "cluster")
                data[[cluster_var]][match(sample_units, data[[ivar]])] else NULL)

    post_periods <- all_periods[all_periods >= g]
    cohort_res[[gkey]] <- list(
      cohort   = g, att = res$att, se = res$se, tstat = res$tstat,
      pvalue   = res$pvalue, ci_lower = res$ci_lower, ci_upper = res$ci_upper,
      n_periods = length(post_periods), n_units = N1, df_resid = res$df_resid
    )
  }

  cohort_df <- as.data.frame(do.call(rbind, lapply(cohort_res, as.data.frame)))
  for (col in names(cohort_df)) cohort_df[[col]] <- as.numeric(cohort_df[[col]])

  # --- Step 3: Overall ATT via pooled regression ---
  # Weights: w_g = N1_g / N_total_treated
  N_treat_total <- sum(cohort_df$n_units)
  w_g <- stats::setNames(cohort_df$n_units / N_treat_total, as.character(cohort_df$cohort))

  # Build pooled Y_bar and D vectors (one obs per unit)
  all_units_all  <- unique(data[[ivar]])
  treated_all    <- unique(data[[ivar]][!never_mask])
  Y_bar_pooled   <- rep(NA_real_, length(all_units_all))
  D_pooled       <- rep(NA_real_, length(all_units_all))
  names(Y_bar_pooled) <- all_units_all
  names(D_pooled)     <- all_units_all

  # Treated units: own-cohort ydot_postavg
  for (g in valid_cohorts) {
    gkey <- as.character(g)
    cohort_units_g <- unique(data[[ivar]][!never_mask & data[[gvar]] == g])
    for (uid in cohort_units_g) {
      if (!is.na(ydot_g_list[[gkey]][as.character(uid)])) {
        Y_bar_pooled[as.character(uid)] <- ydot_g_list[[gkey]][as.character(uid)]
        D_pooled[as.character(uid)]     <- 1
      }
    }
  }

  # Never-treated units: weighted average across cohorts
  for (uid in never_units) {
    uid_c <- as.character(uid)
    w_vals <- numeric(0); y_vals <- numeric(0)
    for (g in valid_cohorts) {
      gkey <- as.character(g)
      val  <- ydot_g_list[[gkey]][uid_c]
      if (!is.null(val) && !is.na(val)) {
        w_vals <- c(w_vals, w_g[gkey])
        y_vals <- c(y_vals, val)
      }
    }
    if (length(y_vals) > 0) {
      wsum <- sum(w_vals)
      Y_bar_pooled[uid_c] <- sum(w_vals * y_vals) / wsum
      D_pooled[uid_c]     <- 0
    }
  }

  # Keep only units with valid Y_bar and D
  valid_mask_pool <- !is.na(Y_bar_pooled) & !is.na(D_pooled)
  Y_reg <- Y_bar_pooled[valid_mask_pool]
  D_reg <- D_pooled[valid_mask_pool]
  unit_ids_reg <- names(Y_reg)

  X_pool <- cbind(intercept = 1, d = D_reg)
  fit_pool <- .ols_fit(X_pool, Y_reg)
  cvec_pool <- if (!is.null(vce) && tolower(vce) == "cluster")
    data[[cluster_var]][match(unit_ids_reg, data[[ivar]])] else NULL
  res_pool <- .extract_inference(fit_pool, vce, cvec_pool)

  att_overall_list <- list(
    att      = res_pool$att, se     = res_pool$se,
    tstat    = res_pool$tstat, pvalue = res_pool$pvalue,
    ci_lower = res_pool$ci_lower, ci_upper = res_pool$ci_upper,
    weights  = w_g, n_cohorts = length(valid_cohorts)
  )

  # --- Step 4: (g, r)-specific effects ---
  gt_results <- list()

  for (g in valid_cohorts) {
    pre_periods  <- all_periods[all_periods < g]
    post_periods <- all_periods[all_periods >= g]
    tpost1       <- min(post_periods)

    if (length(pre_periods) < 1) next
    if (rolling %in% c("detrend","detrendq") && length(pre_periods) < 2) next

    treated_units <- unique(data[[ivar]][!never_mask & data[[gvar]] == g])
    if (control_group == "never_treated") {
      control_units <- never_units
    } else {
      nyt_units     <- unique(data[[ivar]][!never_mask & data[[gvar]] > g])
      control_units <- unique(c(never_units, nyt_units))
    }

    cohort_units <- c(treated_units, control_units)
    df_c <- data[data[[ivar]] %in% cohort_units, , drop = FALSE]
    df_c$post_g <- as.integer(df_c[[tvar]] >= g)
    df_c$d_g    <- as.integer(df_c[[ivar]] %in% treated_units)

    df_t <- apply_transform(df_c, y = y, ivar = ivar, tindex = tvar,
                            post = "post_g", rolling = rolling,
                            tpost1 = tpost1, season_var = season_var)

    per_res <- estimate_period_effects(
      df = df_t, d = "d_g", tindex = tvar,
      post_periods = post_periods, vce = vce,
      cluster_var = cluster_var, controls = controls)

    for (ri in seq_len(nrow(per_res))) {
      r     <- per_res$tindex[ri]
      samp_r <- df_t[df_t[[tvar]] == r, , drop = FALSE]
      gt_results[[length(gt_results) + 1]] <- c(
        list(cohort = g, period = r, event_time = r - g,
             n_treated = sum(samp_r$d_g == 1),
             n_control = sum(samp_r$d_g == 0),
             n_total   = nrow(samp_r)),
        as.list(per_res[ri, c("att","se","tstat","pvalue","ci_lower","ci_upper")])
      )
    }
  }

  gt_df <- if (length(gt_results) > 0) {
    df <- as.data.frame(do.call(rbind, lapply(gt_results, as.data.frame)))
    for (col in names(df)) df[[col]] <- as.numeric(df[[col]])
    df
  } else NULL

  list(
    att_overall        = att_overall_list,
    att_by_cohort      = cohort_df,
    att_by_cohort_time = gt_df,
    cohorts            = valid_cohorts,
    rolling            = rolling,
    control_group      = control_group,
    aggregate          = aggregate,
    vce                = vce
  )
}
