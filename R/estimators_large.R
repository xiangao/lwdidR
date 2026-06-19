# Large-N per-cell estimators for LWDID: RA / IPW / IPWRA.
#
# Faithful port of the point-estimation block in Stata lwdid v2.4.2 `lwdid_large`
# (lines ~1505-1592 of lwdid.ado). For a given cohort g and calendar period t,
# the cross-section is the eligible control sample at time t; the treatment
# indicator dvar_g = 1[gvar == g]. The ATT(g,t) is the coefficient on dvar_g in:
#   ra    : ydot_g ~ dvar_g [+ xc + dvar_g:xc]                 (vce robust)
#   ipw   : ydot_g ~ dvar_g                  [aw = ipw]        (vce robust)
#   ipwra : ydot_g ~ dvar_g [+ xc + dvar_g:xc] [aw = ipw]      (vce robust)
# where xc are the COHORT-CENTERED covariates (v - mean(v among dvar_g==1)), the
# logit propensity score is fit on the RAW covariates, and the ATT weights are
#   ipw = 1 if treated; phat/(1-phat) if control.
# OLS coefficients are invariant to overall weight scaling, so Stata's `aw`
# sum-to-N normalisation does not affect the point estimate (it matters only for
# the influence-function variance handled in Phase 2).

#' Fit one large-N ATT(g,t) cell
#'
#' @param ydot Numeric vector: cohort-g transformed outcome on the cell cross-section.
#' @param dvar Numeric 0/1 treatment indicator (1 = first-treated in cohort g).
#' @param xc Matrix of cohort-centered covariates (cell rows x p), or NULL.
#' @param xraw Matrix of RAW covariates for the propensity score (cell rows x p), or NULL.
#' @param method "ra", "ipw", or "ipwra".
#' @return list(ok, att, ngt, fit, w, used) where `used` is the logical vector of
#'   rows entering estimation. `ok = FALSE` marks an unusable cell (skipped, as in
#'   Stata's `cap ... continue`).
#' @keywords internal
.fit_cell_large <- function(ydot, dvar, xc, xraw, method) {
  bad <- list(ok = FALSE, att = NA_real_, ngt = 0L)

  # propensity score + ATT weights for ipw / ipwra (logit on RAW covariates)
  w <- NULL
  phat <- NULL
  psamp <- rep(TRUE, length(ydot))
  if (method %in% c("ipw", "ipwra")) {
    if (is.null(xraw)) return(bad)  # propensity model needs covariates
    ps_ok <- stats::complete.cases(dvar, xraw)
    fit_l <- tryCatch(
      suppressWarnings(stats::glm(dvar[ps_ok] ~ xraw[ps_ok, , drop = FALSE],
                                  family = stats::binomial())),
      error = function(e) NULL)
    if (is.null(fit_l) || !fit_l$converged) return(bad)
    phat <- rep(NA_real_, length(ydot))
    phat[ps_ok] <- stats::fitted(fit_l)
    # control odds weights; treated weight 1. Guard against phat -> 1 (Inf).
    wt <- ifelse(dvar == 1, 1, phat / (1 - phat))
    if (any(!is.finite(wt[ps_ok]))) return(bad)
    w <- wt
    psamp <- ps_ok
  }

  # The outcome regression includes covariates only for RA and IPWRA. For pure
  # IPW the outcome model is y ~ d (weighted); covariates enter only through the
  # propensity score (matching Stata: reg yvar dvar_g [aw=ipw]).
  xc_out <- if (method == "ipw") NULL else xc

  # assemble estimation sample (listwise deletion, mirroring regress/aw)
  use <- !is.na(ydot) & !is.na(dvar) & psamp
  if (!is.null(xc_out)) use <- use & stats::complete.cases(xc_out)
  if (sum(use) < 2L || length(unique(dvar[use])) < 2L) return(bad)

  y <- ydot[use]; d <- dvar[use]
  wt <- if (is.null(w)) NULL else w[use]

  # design: include centered covariates and their interaction with d (ra/ipwra)
  fit <- tryCatch({
    if (is.null(xc_out)) {
      stats::lm(y ~ d, weights = wt)
    } else {
      X <- xc_out[use, , drop = FALSE]
      stats::lm(y ~ d * X, weights = wt)
    }
  }, error = function(e) NULL)
  if (is.null(fit)) return(bad)

  cf <- stats::coef(fit)
  if (!("d" %in% names(cf)) || is.na(cf[["d"]])) return(bad)

  # stash used-row data for the Phase-2 influence functions
  pU <- if (is.null(phat)) NULL else phat[use]
  dat <- list(
    y    = y,
    d    = d,
    xc   = if (is.null(xc_out)) NULL else xc_out[use, , drop = FALSE],
    xraw = if (is.null(xraw)) NULL else xraw[use, , drop = FALSE],
    w    = wt,
    p    = pU
  )

  list(ok = TRUE,
       att = unname(cf[["d"]]),
       ngt = sum(d == 1L),
       fit = fit, w = wt, used = use, dat = dat)
}
