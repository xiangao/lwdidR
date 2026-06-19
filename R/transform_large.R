# Large-N rolling transformation for LWDID.
#
# Faithful port of the transformation block in Stata lwdid v2.4.2 `lwdid_large`
# (lines ~1356-1465 of lwdid.ado). For each treatment cohort g, each unit is
# residualised by removing its OWN pre-treatment mean (demean) or linear trend
# (detrend), fitted on the unit's pre-window, then applied to ALL periods. The
# transformed outcome is then ANCHORED to 0 at the normalisation period(s):
#   demean : r == -1
#   detrend: r %in% c(-2, -1)        where r = t - g
# The anchoring is mechanical (a baseline), not a fitted value, and it
# propagates into the WATT(r) series, so it must be reproduced exactly.

#' Large-N rolling transformation for one cohort
#'
#' @param df Long panel data frame (one row per unit-period).
#' @param y Outcome column name.
#' @param ivar Unit id column name.
#' @param tvar Calendar time column name (integer-valued).
#' @param g Cohort (first-treatment period) to transform for.
#' @param rolling "demean" or "detrend".
#' @param pre Integer. Number of most-recent pre-periods to use; -1 uses all
#'   pre-periods (t < g). With "detrend", an effective window of >= 2 is required.
#' @return Numeric vector aligned to the rows of `df`: the transformed outcome
#'   `ydot_g` (NA where the unit has no usable pre-window), with anchor cells
#'   set to exactly 0.
#' @keywords internal
.transform_cohort_large <- function(df, y, ivar, tvar, g, rolling, pre = -1L) {
  yv <- df[[y]]
  iv <- df[[ivar]]
  tv <- df[[tvar]]
  n  <- nrow(df)

  # pre-window mask: t < g, optionally restricted to the `pre` most recent
  # pre-periods [g - pre, g - 1].
  in_pre <- tv < g
  if (pre > 0L) in_pre <- in_pre & tv >= (g - pre)
  in_pre <- in_pre & !is.na(yv) & !is.na(tv)

  ydot <- rep(NA_real_, n)

  if (rolling == "demean") {
    # per-unit pre-window mean
    sy <- tapply(yv[in_pre], iv[in_pre], sum)
    ny <- tapply(yv[in_pre], iv[in_pre], length)
    m  <- sy / ny                                  # named by unit
    key <- as.character(iv)
    mu  <- m[key]
    ok  <- !is.na(mu)
    ydot[ok] <- yv[ok] - mu[ok]

  } else if (rolling == "detrend") {
    # per-unit pre-window linear fit y ~ 1 + t, then residual at all t
    pre_df <- data.frame(i = iv[in_pre], t = tv[in_pre], y = yv[in_pre])
    agg <- function(f) tapply(seq_len(nrow(pre_df)), pre_df$i,
                              function(ix) f(pre_df$t[ix], pre_df$y[ix]))
    nn  <- tapply(pre_df$y, pre_df$i, length)
    St  <- tapply(pre_df$t, pre_df$i, sum)
    Stt <- tapply(pre_df$t^2, pre_df$i, sum)
    Sy  <- tapply(pre_df$y, pre_df$i, sum)
    Sty <- agg(function(t, y) sum(t * y))
    denom <- nn * Stt - St^2
    bP <- ifelse(abs(denom) < 1e-12, 0, (nn * Sty - St * Sy) / denom)
    aP <- (Sy - bP * St) / nn
    bP[nn < 2] <- NA_real_; aP[nn < 2] <- NA_real_   # need >= 2 pre-periods
    key <- as.character(iv)
    a <- aP[key]; b <- bP[key]
    ok <- !is.na(a) & !is.na(b)
    ydot[ok] <- yv[ok] - (a[ok] + b[ok] * tv[ok])

  } else {
    stop(sprintf("Large-N rolling must be 'demean' or 'detrend'; got '%s'.", rolling))
  }

  # anchoring to 0 at the normalisation relative period(s)
  r <- tv - g
  anchor_r <- if (rolling == "demean") -1L else c(-2L, -1L)
  ydot[r %in% anchor_r] <- 0

  ydot
}
