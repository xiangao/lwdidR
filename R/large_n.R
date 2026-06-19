# Large-N path driver for LWDID. Faithful port of `lwdid_large` in Stata lwdid
# v2.4.2 (lines ~1174-2466 of lwdid.ado), Phase 1 = point estimates:
# for each treated cohort g and calendar period t, residualise (demean/detrend
# with anchoring), build the eligible control cross-section at t, estimate the
# ATT(g,t) cell via RA/IPW/IPWRA, then aggregate to WATT(r) + Pre/Post averages.
# Influence-function inference (wild cluster bootstrap, sup-t bands) is added in
# Phase 2.

#' Large-N LWDID (point estimates)
#'
#' @param data Long-format panel data frame.
#' @param y Outcome column name.
#' @param ivar Unit id column name.
#' @param tvar Calendar time column name (integer-valued).
#' @param gvar First-treatment-period column (0 or NA = never treated).
#' @param rolling "demean" or "detrend".
#' @param method "ra", "ipw", or "ipwra".
#' @param controls Character vector of covariate names, or NULL.
#' @param pre Integer; number of most-recent pre-periods used in the transform
#'   (-1 = all pre-periods).
#' @param never Logical; if TRUE, controls are never-treated units only.
#' @param level Confidence level (default 95).
#' @return list(attgt = cell-level data frame, watt = aggregated data frame,
#'   cells = list of per-cell fit objects for Phase-2 inference, meta = list).
#' @keywords internal
lwdid_large <- function(data, y, ivar, tvar, gvar, rolling = "demean",
                        method = "ra", controls = NULL, pre = -1L,
                        never = FALSE, level = 95) {
  if (!rolling %in% c("demean", "detrend"))
    stop("Large-N rolling must be 'demean' or 'detrend'.")
  if (!method %in% c("ra", "ipw", "ipwra"))
    stop("Large-N method must be 'ra', 'ipw', or 'ipwra'.")
  if (method %in% c("ipw", "ipwra") && is.null(controls))
    stop("method = '", method, "' requires covariates (controls).")

  # analysis sample: complete cases on the model variables (Stata `touse`)
  core <- c(y, ivar, tvar, gvar, controls)
  touse <- stats::complete.cases(data[, core, drop = FALSE])
  d <- data[touse, , drop = FALSE]

  gv <- d[[gvar]]
  gv[is.na(gv)] <- 0                       # never-treated coded 0
  tv <- d[[tvar]]
  cohorts <- sort(unique(gv[gv > 0]))
  if (length(cohorts) == 0L) stop("No treated cohorts found in gvar.")
  ts <- sort(unique(tv))

  xraw_all <- if (is.null(controls)) NULL else as.matrix(d[, controls, drop = FALSE])

  cells <- list()
  rows  <- list()
  for (g in cohorts) {
    ydot_g <- .transform_cohort_large(d, y, ivar, tvar, g, rolling, pre)

    # cohort-centered covariates: v - mean(v among first-treated-in-g)
    xc_all <- NULL
    if (!is.null(controls)) {
      treated_g <- gv == g
      xc_all <- xraw_all
      for (j in seq_along(controls)) {
        mj <- mean(xraw_all[treated_g, j], na.rm = TRUE)
        xc_all[, j] <- xraw_all[, j] - mj
      }
    }

    for (tt in ts) {
      r <- tt - g
      if (tt < g) {
        cont <- if (never) (gv == 0 | gv == g) else (gv == 0 | gv == g | gv > g)
      } else {
        cont <- if (never) (gv == 0 | gv == g) else (gv == 0 | gv == g | (gv > g & gv > tt))
      }
      cell <- (tv == tt) & cont
      if (!any(cell)) next

      fit <- .fit_cell_large(
        ydot = ydot_g[cell],
        dvar = as.numeric(gv[cell] == g),
        xc   = if (is.null(xc_all)) NULL else xc_all[cell, , drop = FALSE],
        xraw = if (is.null(xraw_all)) NULL else xraw_all[cell, , drop = FALSE],
        method = method)
      if (!isTRUE(fit$ok)) next

      rows[[length(rows) + 1L]] <- data.frame(g = g, t = tt, r = r,
                                              att = fit$att, ngt = fit$ngt)
      # keep cell context for Phase-2 influence functions
      fit$g <- g; fit$t <- tt; fit$r <- r; fit$cell <- which(cell)
      cells[[length(cells) + 1L]] <- fit
    }
  }

  if (length(rows) == 0L) stop("No estimable ATT(g,t) cells.")
  attgt <- do.call(rbind, rows)
  watt  <- .aggregate_large(attgt, rolling)

  list(attgt = attgt, watt = watt, cells = cells, data = d,
       meta = list(rolling = rolling, method = method, controls = controls,
                   pre = pre, never = never, level = level,
                   ivar = ivar, n_units = length(unique(d[[ivar]]))))
}


# User-facing wrapper: run the large-N path, add inference, build a clean object.
#' @keywords internal
.lwdid_large_dispatch <- function(data, y, ivar, tvar, gvar, rolling, method,
                                  controls, pre, never, attgt, ydot,
                                  cluster_var, reps, seed, level, nose) {
  res <- lwdid_large(data, y, ivar, tvar, gvar, rolling = rolling,
                     method = method, controls = controls, pre = pre,
                     never = never, level = level)
  if (!nose) {
    res <- .inference_large(res, res$data, ivar = ivar,
                            cluster_var = cluster_var, reps = reps,
                            seed = seed, level = level)
  }

  # optional transformed outcomes (one y{g}d column per cohort)
  ydot_df <- NULL
  if (isTRUE(ydot)) {
    d <- res$data
    gv <- d[[gvar]]; gv[is.na(gv)] <- 0
    cohorts <- sort(unique(gv[gv > 0]))
    ydot_df <- d[, c(ivar, tvar, gvar), drop = FALSE]
    for (g in cohorts) {
      ydot_df[[paste0("y", g, "d")]] <-
        .transform_cohort_large(d, y, ivar, tvar, g, rolling, pre)
    }
  }

  structure(
    list(design = "large_n",
         watt = res$watt,
         attgt = if (isTRUE(attgt)) res$attgt else NULL,
         ydot = ydot_df,
         rolling = rolling, method = method, controls = controls,
         pre = pre, never = never, level = level, nose = nose,
         n_units = res$meta$n_units,
         y = y, ivar = ivar, tvar = tvar, gvar = gvar),
    class = c("lwdid_large", "lwdid"))
}
