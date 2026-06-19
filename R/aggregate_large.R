# Aggregate large-N ATT(g,t) cells to the event-study WATT(r) path and the
# Pre/Post averages. Faithful port of Stage 2b in Stata lwdid v2.4.2
# (lines ~1866-1936 of lwdid.ado).
#
# WATT(r) is the TREATED-COUNT-weighted average of ATT(g,t) over cells sharing
# the same relative period r = t - g:  weight = Ngt / sum(Ngt).
# Anchoring (demean: r = -1; detrend: r in {-2,-1}) forces WATT(r) = 0.
# Pre_avg averages pre cells (r < -1 demean, r < -2 detrend); Post_avg averages
# post cells (r >= 0); weights are Ngt normalised within each block.

#' Aggregate ATT(g,t) cells to WATT(r) + Pre/Post averages
#'
#' @param attgt Data frame with columns g, t, r, att, ngt (one row per cell).
#' @param rolling "demean" or "detrend" (controls the anchor periods).
#' @return Data frame with columns effect, ryear, watt, n_cells, n_units.
#'   Rows: Pre_avg, Post_avg, then one row per relative period r (sorted).
#'   SE / t / p / CI columns are added later by the inference layer (Phase 2).
#' @keywords internal
.aggregate_large <- function(attgt, rolling) {
  anchor_r <- if (rolling == "demean") -1L else c(-2L, -1L)
  pre_max  <- if (rolling == "demean") -2L else -3L   # pre block: r <= pre_max

  wsum <- function(att, ngt) {
    keep <- !is.na(att) & !is.na(ngt) & ngt > 0
    if (!any(keep)) return(c(val = NA_real_, n_cells = 0, n_units = 0))
    a <- att[keep]; n <- ngt[keep]
    c(val = sum((n / sum(n)) * a), n_cells = sum(keep), n_units = sum(n))
  }

  # --- WATT(r) for each relative period present ---
  rs <- sort(unique(attgt$r))
  watt_rows <- lapply(rs, function(rr) {
    cells <- attgt[attgt$r == rr, , drop = FALSE]
    if (rr %in% anchor_r) {
      data.frame(effect = sprintf("WATT(%d)", rr), ryear = rr, watt = 0,
                 n_cells = nrow(cells), n_units = sum(cells$ngt, na.rm = TRUE))
    } else {
      s <- wsum(cells$att, cells$ngt)
      data.frame(effect = sprintf("WATT(%d)", rr), ryear = rr, watt = unname(s["val"]),
                 n_cells = unname(s["n_cells"]), n_units = unname(s["n_units"]))
    }
  })
  watt <- do.call(rbind, watt_rows)

  # --- Pre / Post averages ---
  pre_cells  <- attgt[attgt$r <= pre_max, , drop = FALSE]
  post_cells <- attgt[attgt$r >= 0,       , drop = FALSE]
  sp <- wsum(pre_cells$att,  pre_cells$ngt)
  so <- wsum(post_cells$att, post_cells$ngt)
  avg <- data.frame(
    effect  = c("Pre_avg", "Post_avg"),
    ryear   = c(NA_real_, NA_real_),
    watt    = c(unname(sp["val"]), unname(so["val"])),
    n_cells = c(unname(sp["n_cells"]), unname(so["n_cells"])),
    n_units = c(unname(sp["n_units"]), unname(so["n_units"]))
  )

  out <- rbind(avg, watt)
  rownames(out) <- NULL
  out
}
