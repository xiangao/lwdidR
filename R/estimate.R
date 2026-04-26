# Cross-sectional OLS estimation for LWDID
# Mirrors Python lwdid/staggered/estimation.py and aggregation.py logic


#' Prepare centred controls and D*X interactions
#'
#' @param df_sample Cross-section data frame
#' @param d Treatment indicator column name
#' @param controls Character vector of control variable names (or NULL)
#' @return List with: include, X_centered, interactions, RHS_varnames
#' @keywords internal
prepare_controls <- function(df_sample, d, controls) {
  if (is.null(controls) || length(controls) == 0) {
    return(list(include = FALSE, X_centered = NULL, interactions = NULL,
                RHS_varnames = character(0)))
  }
  K  <- length(controls)
  N1 <- sum(df_sample[[d]] == 1, na.rm = TRUE)
  N0 <- sum(df_sample[[d]] == 0, na.rm = TRUE)
  if (N1 > K + 1 && N0 > K + 1) {
    xbar <- vapply(controls, function(x)
      mean(df_sample[[x]][df_sample[[d]] == 1], na.rm = TRUE), numeric(1))
    X_centered <- as.data.frame(
      lapply(controls, function(x) df_sample[[x]] - xbar[x]))
    names(X_centered) <- paste0(controls, "_c")
    interactions <- as.data.frame(
      lapply(controls, function(x)
        df_sample[[d]] * X_centered[[paste0(x, "_c")]]))
    names(interactions) <- paste0("d_", controls, "_c")
    list(include = TRUE, X_centered = X_centered,
         interactions = interactions, xbar = xbar,
         RHS_varnames = c(controls, paste0("d_", controls, "_c")))
  } else {
    warning(sprintf(
      "Controls not applied: N1=%d, N0=%d, K+1=%d. Controls ignored.", N1, N0, K + 1))
    list(include = FALSE, X_centered = NULL, interactions = NULL,
         RHS_varnames = character(0))
  }
}


# Build design matrix from data frame, treatment column, and optional controls
.build_X <- function(samp, d, controls, ctrl_spec) {
  if (!is.null(ctrl_spec) && ctrl_spec$include) {
    cbind(intercept = 1,
          d = as.numeric(samp[[d]]),
          as.matrix(samp[, controls, drop = FALSE]),
          as.matrix(ctrl_spec$interactions))
  } else {
    cbind(intercept = 1, d = as.numeric(samp[[d]]))
  }
}


# OLS fit from design matrix X and response y
# Returns list: coef, residuals, fitted, df_resid, X, y
.ols_fit <- function(X, y) {
  b    <- as.numeric(solve(t(X) %*% X, t(X) %*% y))
  yhat <- as.numeric(X %*% b)
  e    <- y - yhat
  list(coef = b, residuals = e, fitted = yhat,
       df_resid = nrow(X) - ncol(X), X = X, y = y)
}


# Compute SE vector for OLS fit given vce type
# Returns named numeric vector (one SE per coefficient)
.compute_se <- function(fit, vce, cluster_vec = NULL) {
  X    <- fit$X
  e    <- fit$residuals
  n    <- nrow(X)
  k    <- ncol(X)
  XtX_inv <- solve(t(X) %*% X)

  if (is.null(vce)) {
    # Homoskedastic OLS
    sigma2 <- sum(e^2) / (n - k)
    vcov_mat <- sigma2 * XtX_inv
  } else if (tolower(vce) %in% c("hc1", "robust")) {
    # HC1: n/(n-k) correction
    meat <- t(X) %*% diag(e^2) %*% X
    vcov_mat <- (n / (n - k)) * XtX_inv %*% meat %*% XtX_inv
  } else if (tolower(vce) == "hc3") {
    # HC3: residuals divided by (1 - h_ii)^2
    # h_ii = X_i (X'X)^{-1} X_i'
    h  <- rowSums((X %*% XtX_inv) * X)   # leverage values
    # Clamp h to avoid divide by zero (treat 0/0 as 0)
    e_adj <- ifelse(abs(1 - h) < 1e-10, 0, e / (1 - h))
    meat <- t(X) %*% diag(e_adj^2) %*% X
    vcov_mat <- XtX_inv %*% meat %*% XtX_inv
  } else if (tolower(vce) == "cluster") {
    if (is.null(cluster_vec)) stop("vce='cluster' requires cluster_vec")
    g_ids <- unique(cluster_vec)
    G     <- length(g_ids)
    meat  <- matrix(0, k, k)
    for (g in g_ids) {
      idx <- cluster_vec == g
      sc  <- colSums(X[idx, , drop = FALSE] * e[idx])
      meat <- meat + outer(sc, sc)
    }
    correction <- (G / (G - 1)) * ((n - 1) / (n - k))
    vcov_mat <- correction * XtX_inv %*% meat %*% XtX_inv
  } else {
    stop(sprintf(
      "Unknown vce='%s'. Use NULL, 'hc1', 'hc3', 'cluster', 'wildboot', 'permutation'.",
      vce))
  }
  sqrt(pmax(0, diag(vcov_mat)))
}


# Wild cluster bootstrap inference
# Uses Rademacher weights (one per cluster or per unit).
# Returns same list structure as .extract_inference.
.wild_bootstrap <- function(fit, vce_inner = "hc3", nboot = 999,
                             cluster_vec = NULL, alpha = 0.05) {
  X       <- fit$X
  yhat    <- fit$fitted
  e       <- fit$residuals
  n       <- nrow(X)
  att_obs <- fit$coef[2]

  se_obs <- .compute_se(fit, vce_inner, cluster_vec)[2]
  t_obs  <- att_obs / se_obs

  groups     <- if (!is.null(cluster_vec)) cluster_vec else seq_len(n)
  unique_grps <- unique(groups)
  G           <- length(unique_grps)

  b_star <- numeric(nboot)
  t_star <- numeric(nboot)
  for (s in seq_len(nboot)) {
    w_grp        <- sample(c(-1, 1), G, replace = TRUE)
    names(w_grp) <- as.character(unique_grps)
    w            <- w_grp[as.character(groups)]
    fit_star     <- .ols_fit(X, yhat + w * e)
    se_star      <- .compute_se(fit_star, vce_inner, cluster_vec)[2]
    b_star[s]    <- fit_star$coef[2]
    t_star[s]    <- (fit_star$coef[2] - att_obs) / se_star
  }

  list(
    att      = unname(att_obs),
    se       = stats::sd(b_star),
    tstat    = unname(t_obs),
    pvalue   = mean(abs(t_star) >= abs(t_obs)),
    ci_lower = unname(stats::quantile(b_star, alpha / 2)),
    ci_upper = unname(stats::quantile(b_star, 1 - alpha / 2)),
    df_resid = fit$df_resid,
    N        = n
  )
}


# Permutation (randomisation) inference
# Permutes the treatment column; rebuilds D*X interactions when controls present.
.permutation_test <- function(X, y, d_col_idx = 2L, nperm = 999, alpha = 0.05) {
  n       <- nrow(X)
  k       <- ncol(X)
  fit_obs <- .ols_fit(X, y)
  att_obs <- fit_obs$coef[d_col_idx]
  K_ctrl  <- (k - 2L) %/% 2L   # number of controls (0 if none)

  att_star <- numeric(nperm)
  for (s in seq_len(nperm)) {
    X_perm            <- X
    X_perm[, d_col_idx] <- X[sample.int(n), d_col_idx]
    if (K_ctrl > 0L) {
      ctrl_start  <- 3L
      inter_start <- 3L + K_ctrl
      for (ki in seq_len(K_ctrl)) {
        X_perm[, inter_start + ki - 1L] <-
          X_perm[, d_col_idx] * X[, ctrl_start + ki - 1L]
      }
    }
    att_star[s] <- .ols_fit(X_perm, y)$coef[d_col_idx]
  }

  se_perm <- stats::sd(att_star)
  z_crit  <- stats::qnorm(1 - alpha / 2)
  list(
    att      = unname(att_obs),
    se       = se_perm,
    tstat    = unname(att_obs / se_perm),
    pvalue   = mean(abs(att_star) >= abs(att_obs)),
    ci_lower = unname(att_obs - z_crit * se_perm),
    ci_upper = unname(att_obs + z_crit * se_perm),
    df_resid = fit_obs$df_resid,
    N        = n
  )
}


# Extract inference: att (coef[2]), se, tstat, pvalue, CI
.extract_inference <- function(fit, vce, cluster_vec = NULL, alpha = 0.05,
                                nboot = 999, nperm = 999, vce_inner = "hc3") {
  if (!is.null(vce) && tolower(vce) == "wildboot")
    return(.wild_bootstrap(fit, vce_inner, nboot, cluster_vec, alpha))
  if (!is.null(vce) && tolower(vce) == "permutation")
    return(.permutation_test(fit$X, fit$y, 2L, nperm, alpha))
  ses <- .compute_se(fit, vce, cluster_vec)

  if (!is.null(vce) && tolower(vce) == "cluster") {
    G  <- length(unique(cluster_vec))
    df <- G - 1
  } else {
    df <- fit$df_resid
  }

  att   <- fit$coef[2]
  se    <- ses[2]
  tstat <- att / se
  pval  <- 2 * stats::pt(-abs(tstat), df = df)
  tcrit <- stats::qt(1 - alpha / 2, df = df)
  list(att = unname(att), se = unname(se), tstat = unname(tstat),
       pvalue = unname(pval),
       ci_lower = unname(att - tcrit * se),
       ci_upper = unname(att + tcrit * se),
       df_resid = df, N = nrow(fit$X))
}


#' Estimate ATT via OLS on firstpost cross-section (uses ydot_postavg)
#'
#' @param df Transformed data with ydot_postavg and firstpost columns
#' @param d Treatment indicator column
#' @param vce Variance estimator
#' @param cluster_var Cluster variable name
#' @param controls Control variable names
#' @return Named list with att, se, tstat, pvalue, ci_lower, ci_upper, N, df_resid
#' @keywords internal
estimate_att <- function(df, d, vce = NULL, cluster_var = NULL, controls = NULL,
                         nboot = 999, nperm = 999, vce_inner = "hc3") {
  samp <- df[df$firstpost == TRUE, , drop = FALSE]
  N    <- nrow(samp)
  N1   <- sum(samp[[d]] == 1, na.rm = TRUE)
  N0   <- N - N1

  if (N < 3) stop("Insufficient sample (N < 3).")
  if (N1 == 0) stop("No treated units in firstpost sample.")
  if (N0 == 0) stop("No control units in firstpost sample.")

  ctrl_spec <- prepare_controls(samp, d, controls)
  X    <- .build_X(samp, d, controls, ctrl_spec)
  y    <- samp$ydot_postavg
  fit  <- .ols_fit(X, y)

  cvec <- if (!is.null(vce) && tolower(vce) == "cluster") samp[[cluster_var]] else NULL
  .extract_inference(fit, vce, cvec, nboot = nboot, nperm = nperm, vce_inner = vce_inner)
}


#' Estimate period-specific ATTs (uses ydot)
#'
#' @param df Transformed panel data
#' @param d Treatment indicator column
#' @param tindex Time index column
#' @param post_periods Integer vector of post-period indices
#' @param vce Variance estimator
#' @param cluster_var Cluster variable name
#' @param controls Control variable names
#' @return Data frame with: tindex, att, se, tstat, pvalue, ci_lower, ci_upper, N
#' @keywords internal
estimate_period_effects <- function(df, d, tindex, post_periods,
                                    vce = NULL, cluster_var = NULL,
                                    controls = NULL,
                                    nboot = 999, nperm = 999, vce_inner = "hc3") {
  rows <- lapply(post_periods, function(t) {
    samp <- df[df[[tindex]] == t, , drop = FALSE]
    N    <- nrow(samp)
    N1   <- sum(samp[[d]] == 1, na.rm = TRUE)

    empty_row <- data.frame(tindex = t, att = NA_real_, se = NA_real_,
                            tstat = NA_real_, pvalue = NA_real_,
                            ci_lower = NA_real_, ci_upper = NA_real_, N = N)

    if (N < 3 || N1 == 0 || (N - N1) == 0) return(empty_row)

    tryCatch({
      ctrl_spec <- prepare_controls(samp, d, controls)
      X    <- .build_X(samp, d, controls, ctrl_spec)
      y    <- samp$ydot
      fit  <- .ols_fit(X, y)
      cvec <- if (!is.null(vce) && tolower(vce) == "cluster") samp[[cluster_var]] else NULL
      res  <- .extract_inference(fit, vce, cvec,
                                  nboot = nboot, nperm = nperm, vce_inner = vce_inner)
      data.frame(tindex = t, att = res$att, se = res$se,
                 tstat = res$tstat, pvalue = res$pvalue,
                 ci_lower = res$ci_lower, ci_upper = res$ci_upper, N = N)
    }, error = function(e) {
      warning(sprintf("Period %d regression failed: %s", t, conditionMessage(e)))
      empty_row
    })
  })
  do.call(rbind, rows)
}
