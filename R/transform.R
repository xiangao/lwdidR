# Unit-specific pre-treatment transformations for LWDID
# Mirrors Python lwdid/transformations.py logic


#' Apply unit-specific transformation to panel data
#'
#' Computes ydot (residualised outcome), ydot_postavg (post-treatment average
#' of ydot per unit), and marks the firstpost cross-section.
#'
#' @param df Data frame (long format, one row per unit-period)
#' @param y Outcome column name
#' @param ivar Unit identifier column name
#' @param tindex Integer time index column name
#' @param post Binary post-treatment indicator column name (0=pre, 1=post)
#' @param rolling Transformation method: "demean", "detrend", "demeanq", "detrendq"
#' @param tpost1 First post-treatment period index
#' @param season_var Column name of seasonal indicator (required for demeanq/detrendq)
#' @return Data frame with added columns: ydot, ydot_postavg, firstpost
#' @export
apply_transform <- function(df, y, ivar, tindex, post, rolling, tpost1,
                            season_var = NULL) {
  df <- df  # work on copy (R passes by value)

  units <- unique(df[[ivar]])
  df$ydot <- NA_real_

  if (rolling == "demean") {
    for (uid in units) {
      idx_pre <- df[[ivar]] == uid & df[[post]] == 0
      idx_all <- df[[ivar]] == uid
      if (sum(idx_pre) < 1) {
        stop(sprintf("Unit %s has no pre-treatment observations (demean requires >= 1).", uid))
      }
      y_pre_mean <- mean(df[[y]][idx_pre], na.rm = TRUE)
      df$ydot[idx_all] <- df[[y]][idx_all] - y_pre_mean
    }

  } else if (rolling == "detrend") {
    for (uid in units) {
      idx_pre <- df[[ivar]] == uid & df[[post]] == 0
      idx_all <- df[[ivar]] == uid
      n_pre <- sum(idx_pre)
      if (n_pre < 2) {
        stop(sprintf("Unit %s has only %d pre-period(s) (detrend requires >= 2).", uid, n_pre))
      }
      y_pre  <- df[[y]][idx_pre]
      t_pre  <- df[[tindex]][idx_pre]
      t_mean <- mean(t_pre)
      fit    <- .lm_simple(y_pre, t_pre - t_mean)
      t_all  <- df[[tindex]][idx_all]
      yhat   <- fit$alpha + fit$beta * (t_all - t_mean)
      df$ydot[idx_all] <- df[[y]][idx_all] - yhat
    }

  } else if (rolling == "demeanq") {
    if (is.null(season_var)) stop("rolling='demeanq' requires season_var.")
    for (uid in units) {
      idx_pre <- df[[ivar]] == uid & df[[post]] == 0
      idx_all <- df[[ivar]] == uid
      unit_data_pre <- df[idx_pre, , drop = FALSE]
      unit_data_all <- df[idx_all, , drop = FALSE]
      res <- .demeanq_unit(unit_data_pre, unit_data_all, y, season_var)
      df$ydot[idx_all] <- res
    }

  } else if (rolling == "detrendq") {
    if (is.null(season_var)) stop("rolling='detrendq' requires season_var.")
    for (uid in units) {
      idx_pre <- df[[ivar]] == uid & df[[post]] == 0
      idx_all <- df[[ivar]] == uid
      unit_data_pre <- df[idx_pre, , drop = FALSE]
      unit_data_all <- df[idx_all, , drop = FALSE]
      res <- .detrendq_unit(unit_data_pre, unit_data_all, y, tindex, season_var)
      df$ydot[idx_all] <- res
    }

  } else {
    stop(sprintf("Invalid rolling method: '%s'. Choose demean, detrend, demeanq, detrendq.", rolling))
  }

  # Post-treatment average of ydot per unit
  post_mask <- df[[post]] == 1
  ydot_postavg <- tapply(df$ydot[post_mask], df[[ivar]][post_mask], mean, na.rm = TRUE)
  df$ydot_postavg <- ydot_postavg[as.character(df[[ivar]])]

  # Mark first post-treatment cross-section
  df$firstpost <- df[[tindex]] == tpost1 & !is.na(df$ydot_postavg)

  df
}


# OLS slope + intercept for simple linear regression (y ~ 1 + t_centered)
.lm_simple <- function(y, t_c) {
  n <- length(y)
  valid <- !is.na(y) & !is.na(t_c)
  y <- y[valid]; t_c <- t_c[valid]
  n <- length(y)
  if (n < 2) return(list(alpha = NA_real_, beta = NA_real_))
  Stt <- sum(t_c^2) - sum(t_c)^2 / n
  Sty <- sum(t_c * y) - sum(t_c) * sum(y) / n
  if (abs(Stt) < 1e-10) return(list(alpha = mean(y), beta = 0))
  beta  <- Sty / Stt
  alpha <- mean(y) - beta * mean(t_c)
  list(alpha = alpha, beta = beta)
}


# Seasonal demeaning for a single unit (demeanq)
# unit_data_pre: pre-period rows; unit_data_all: all rows
.demeanq_unit <- function(unit_data_pre, unit_data_all, y, season_var) {
  y_pre  <- unit_data_pre[[y]]
  s_pre  <- as.factor(unit_data_pre[[season_var]])
  valid  <- !is.na(y_pre) & !is.na(unit_data_pre[[season_var]])
  if (sum(valid) <= nlevels(s_pre)) {
    warning("demeanq: insufficient pre-period obs for reliable estimation. Returning NA.")
    return(rep(NA_real_, nrow(unit_data_all)))
  }
  fit <- stats::lm(y_pre ~ s_pre, na.action = stats::na.omit)
  # Predict for all periods using same factor levels
  s_all  <- factor(unit_data_all[[season_var]], levels = levels(s_pre))
  yhat   <- stats::predict(fit, newdata = data.frame(s_pre = s_all))
  unit_data_all[[y]] - yhat
}


# Seasonal detrending for a single unit (detrendq)
.detrendq_unit <- function(unit_data_pre, unit_data_all, y, tindex, season_var) {
  y_pre  <- unit_data_pre[[y]]
  t_pre  <- unit_data_pre[[tindex]]
  s_pre  <- as.factor(unit_data_pre[[season_var]])
  valid  <- !is.na(y_pre) & !is.na(t_pre) & !is.na(unit_data_pre[[season_var]])
  n_seas <- nlevels(s_pre)
  if (sum(valid) <= (n_seas + 1)) {
    warning("detrendq: insufficient pre-period obs for reliable estimation. Returning NA.")
    return(rep(NA_real_, nrow(unit_data_all)))
  }
  t_mean <- mean(t_pre[valid])
  t_c_pre <- t_pre - t_mean
  fit <- stats::lm(y_pre ~ t_c_pre + s_pre, na.action = stats::na.omit)
  t_c_all <- unit_data_all[[tindex]] - t_mean
  s_all   <- factor(unit_data_all[[season_var]], levels = levels(s_pre))
  yhat    <- stats::predict(fit, newdata = data.frame(t_c_pre = t_c_all, s_pre = s_all))
  unit_data_all[[y]] - yhat
}
