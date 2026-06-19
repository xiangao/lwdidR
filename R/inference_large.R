# Large-N inference for LWDID: influence functions + wild cluster bootstrap +
# simultaneous (sup-t) bands. Faithful port of the Mata inference code in Stata
# lwdid v2.4.2 (per-cell IFs at lines ~1638-1804; aggregation/bootstrap/bands at
# lines ~1985-2188).
#
# Each cell contributes a length-N influence function for its ATT(g,t). These
# are aggregated (treated-count weights) to event-time r and to Pre/Post blocks,
# centered, then a wild *cluster* multiplier (Rademacher) bootstrap on the
# aggregated IFs yields pointwise SEs and a sup-t simultaneous band for WATT(r).
# Pre/Post averages use pointwise normal CIs.
#
# NOTE on validation: bootstrap SEs depend on the RNG, so they cannot match Stata
# bit-for-bit. For Rademacher weights the bootstrap variance converges to the
# analytic cluster sum sum_cl (sum_{i in cl} IF_i)^2, which we ALSO return
# (`se_analytic`) as a deterministic cross-check against Stata's reported SE.

#' Per-cell influence function for the ATT(g,t) coefficient
#'
#' @param cell A fitted-cell object from `.fit_cell_large` (with `$dat`, `$fit`).
#' @param method "ra", "ipw", or "ipwra".
#' @return Numeric vector of length = #used rows of the cell (the IF on the
#'   estimation sample).
#' @keywords internal
.cell_if <- function(cell, method) {
  dat <- cell$dat
  d <- dat$d; y <- dat$y
  n <- length(d)

  if (method == "ra") {
    Z  <- stats::model.matrix(cell$fit)
    u  <- stats::residuals(cell$fit)
    jd <- match("d", colnames(Z))
    ZZinv <- solve(crossprod(Z))
    IFb <- (Z %*% ZZinv) * u            # row i scaled by u_i
    return(IFb[, jd])

  } else if (method == "ipw") {
    att <- cell$att
    meanD <- mean(d)
    p <- dat$p
    plug <- (d / meanD) * (y - att) -
            ((1 - d) * p / ((1 - p) * meanD)) * y
    X <- cbind(1, dat$xraw)
    W <- p * (1 - p)
    Ainv <- solve(crossprod(X, X * W))
    s   <- X * (d - p)                  # n x k
    IFg <- s %*% Ainv
    gvec <- (1 - d) * y * p / (1 - p)
    Gamma <- colMeans(X * gvec)         # k-vector
    corr <- as.numeric(IFg %*% Gamma)
    return((plug - corr) / n)

  } else if (method == "ipwra") {
    Z <- stats::model.matrix(cell$fit)
    u <- stats::residuals(cell$fit)
    w <- dat$w
    jd <- match("d", colnames(Z))
    Qwinv <- solve(crossprod(Z, Z * w))
    M <- Z * (w * u)
    if (is.null(dat$xc)) {
      IFb <- M %*% Qwinv
      return(IFb[, jd])
    }
    p <- dat$p
    X <- cbind(1, dat$xraw)
    Ainv <- solve(crossprod(X, X * (p * (1 - p))))
    S   <- X * (d - p)
    IFg <- S %*% Ainv
    H   <- crossprod(Z * ((1 - d) * w * u), X)   # ncolZ x ncolX
    IFb <- (M + IFg %*% t(H)) %*% Qwinv
    return(IFb[, jd])
  }
  stop("unknown method")
}


#' Add wild-cluster-bootstrap SEs + sup-t bands to a large-N result
#'
#' @param res Output of `lwdid_large`.
#' @param d The analysis-sample data frame (touse subset) used inside `lwdid_large`.
#' @param ivar Unit id column name (default cluster).
#' @param cluster_var Optional clustering column name (defaults to `ivar`).
#' @param reps Bootstrap replications (default 999).
#' @param seed Optional integer seed.
#' @param level Confidence level (default from res$meta).
#' @return `res` with the `watt` data frame extended by se, se_analytic, t,
#'   p_value, lower_ci, upper_ci.
#' @keywords internal
.inference_large <- function(res, d, ivar, cluster_var = NULL,
                             reps = 999L, seed = NULL, level = NULL) {
  rolling <- res$meta$rolling
  method  <- res$meta$method
  if (is.null(level)) level <- res$meta$level
  alpha <- (100 - level) / 100
  zpt   <- stats::qnorm(1 - alpha / 2)

  Nobs <- nrow(d)
  cells <- res$cells
  K <- length(cells)

  # ---- build IF_mat (Nobs x K): each cell's IF scattered to global rows ----
  IF_mat <- matrix(0, Nobs, K)
  cell_g <- numeric(K); cell_t <- numeric(K); cell_r <- numeric(K); cell_n <- numeric(K)
  for (k in seq_len(K)) {
    ck <- cells[[k]]
    if_used <- .cell_if(ck, method)
    grows <- ck$cell[ck$used]            # global rows of the estimation sample
    IF_mat[grows, k] <- if_used
    cell_g[k] <- ck$g; cell_t[k] <- ck$t; cell_r[k] <- ck$r; cell_n[k] <- ck$ngt
  }

  anchor_r <- if (rolling == "demean") -1L else c(-2L, -1L)
  pre_max  <- if (rolling == "demean") -2L else -3L

  # ---- aggregation weights (must mirror .aggregate_large) ----
  rs <- sort(unique(cell_r))
  # IF_r: Nobs x length(rs)
  IF_r <- matrix(0, Nobs, length(rs))
  for (j in seq_along(rs)) {
    rr <- rs[j]
    if (rr %in% anchor_r) next                  # anchored -> 0 column
    sel <- which(cell_r == rr & cell_n > 0)
    if (!length(sel)) next
    wk <- cell_n[sel] / sum(cell_n[sel])
    IF_r[, j] <- IF_mat[, sel, drop = FALSE] %*% wk
  }
  # IF_avg: Nobs x 2 (Pre, Post)
  pre_sel  <- which(cell_r <= pre_max & cell_n > 0)
  post_sel <- which(cell_r >= 0       & cell_n > 0)
  IF_avg <- matrix(0, Nobs, 2)
  if (length(pre_sel))  IF_avg[, 1] <- IF_mat[, pre_sel,  drop = FALSE] %*% (cell_n[pre_sel]  / sum(cell_n[pre_sel]))
  if (length(post_sel)) IF_avg[, 2] <- IF_mat[, post_sel, drop = FALSE] %*% (cell_n[post_sel] / sum(cell_n[post_sel]))

  # ---- center each aggregated IF column ----
  IF_r   <- sweep(IF_r,   2, colMeans(IF_r))
  IF_avg <- sweep(IF_avg, 2, colMeans(IF_avg))

  # ---- cluster structure (default = unit id) ----
  clv <- if (is.null(cluster_var)) d[[ivar]] else d[[cluster_var]]
  cl_idx <- as.integer(factor(clv))
  n_cl <- max(cl_idx)

  # analytic cluster-robust SE (deterministic): sqrt(sum_cl (sum_i IF)^2)
  cluster_se <- function(IF) {
    apply(IF, 2, function(col) {
      cs <- tapply(col, cl_idx, sum)
      sqrt(sum(cs^2))
    })
  }
  se_r_analytic   <- cluster_se(IF_r)
  se_avg_analytic <- cluster_se(IF_avg)

  # ---- wild cluster bootstrap ----
  if (!is.null(seed)) set.seed(seed)
  BS_r   <- matrix(NA_real_, reps, length(rs))
  BS_avg <- matrix(NA_real_, reps, 2)
  for (b in seq_len(reps)) {
    xi_cl <- sample(c(-1, 1), n_cl, replace = TRUE)
    xi_i  <- xi_cl[cl_idx]
    BS_r[b, ]   <- colSums(IF_r   * xi_i)
    BS_avg[b, ] <- colSums(IF_avg * xi_i)
  }
  se_r   <- apply(BS_r,   2, stats::var); se_r   <- sqrt(se_r)
  se_avg <- sqrt(apply(BS_avg, 2, stats::var))

  # sup-t simultaneous critical value over event times
  se_vec <- se_r
  Tstar <- apply(BS_r, 1, function(row) {
    z <- ifelse(se_vec > 0 & is.finite(se_vec), abs(row / se_vec), 0)
    max(z)
  })
  Tsort <- sort(Tstar)
  q_idx <- min(length(Tsort), max(1L, ceiling((1 - alpha) * length(Tsort))))
  c_sup <- Tsort[q_idx]

  # ---- assemble into res$watt ----
  w <- res$watt
  w$se <- NA_real_; w$se_analytic <- NA_real_
  w$t_stat <- NA_real_; w$p_value <- NA_real_
  w$lower_ci <- NA_real_; w$upper_ci <- NA_real_

  for (j in seq_along(rs)) {
    rr <- rs[j]
    row <- which(w$effect == sprintf("WATT(%d)", rr))
    if (!length(row)) next
    if (rr %in% anchor_r) {
      w$se[row] <- 0; w$se_analytic[row] <- 0
      w$t_stat[row] <- NA; w$p_value[row] <- NA
      w$lower_ci[row] <- 0; w$upper_ci[row] <- 0
    } else {
      w$se[row] <- se_r[j]; w$se_analytic[row] <- se_r_analytic[j]
      th <- w$watt[row]
      w$t_stat[row] <- if (se_r[j] > 0) th / se_r[j] else NA
      w$p_value[row] <- if (se_r[j] > 0) 2 * stats::pnorm(-abs(th / se_r[j])) else NA
      w$lower_ci[row] <- th - c_sup * se_r[j]      # simultaneous band
      w$upper_ci[row] <- th + c_sup * se_r[j]
    }
  }
  for (a in 1:2) {
    lab <- c("Pre_avg", "Post_avg")[a]
    row <- which(w$effect == lab)
    if (!length(row)) next
    th <- w$watt[row]
    w$se[row] <- se_avg[a]; w$se_analytic[row] <- se_avg_analytic[a]
    w$t_stat[row] <- if (se_avg[a] > 0) th / se_avg[a] else NA
    w$p_value[row] <- if (se_avg[a] > 0) 2 * stats::pnorm(-abs(th / se_avg[a])) else NA
    w$lower_ci[row] <- th - zpt * se_avg[a]        # pointwise normal
    w$upper_ci[row] <- th + zpt * se_avg[a]
  }

  res$watt <- w
  res$c_sup <- c_sup
  res
}
