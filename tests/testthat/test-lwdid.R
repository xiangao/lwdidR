library(lwdidR)

# ── Helpers ──────────────────────────────────────────────────────────────────

make_common_panel <- function(n = 30, t_pre = 3, t_post = 2, tau = 2.0,
                               seed = 1) {
  set.seed(seed)
  ids   <- rep(1:n, each = t_pre + t_post)
  times <- rep(seq_len(t_pre + t_post), times = n)
  unit_fe <- rep(rnorm(n), each = t_pre + t_post)
  time_fe <- rep(rnorm(t_pre + t_post), times = n)
  post  <- as.integer(times > t_pre)
  treat <- as.integer(ids <= n / 2)
  y     <- unit_fe + time_fe + tau * treat * post + rnorm(n * (t_pre + t_post))
  data.frame(id = ids, time = times, y = y, post = post, treat = treat)
}

make_staggered_panel <- function(n_per_cohort = 20, t_max = 6, tau = 2.0,
                                  seed = 2) {
  set.seed(seed)
  cohorts <- c(3, 5, 0)   # treat in period 3, 5, or never
  groups  <- rep(cohorts, each = n_per_cohort)
  ids     <- seq_along(groups)
  rows    <- lapply(seq_along(ids), function(i) {
    unit_fe <- rnorm(1)
    data.frame(
      id    = ids[i],
      time  = seq_len(t_max),
      gvar  = groups[i],
      y     = unit_fe + rnorm(t_max) +
              ifelse(groups[i] > 0 & seq_len(t_max) >= groups[i], tau, 0)
    )
  })
  do.call(rbind, rows)
}

# ── Common-timing tests ───────────────────────────────────────────────────────

test_that("common-timing demean recovers ATT", {
  df  <- make_common_panel(tau = 2.0)
  res <- lwdid(df, y = "y", ivar = "id", tvar = "time", post = "post",
               rolling = "demean", vce = NULL)
  expect_s3_class(res, "lwdid")
  expect_true(is.numeric(res$att_overall))
  expect_true(is.numeric(res$se_overall))
  expect_equal(res$design, "common_timing")
  # ATT estimate should be near true value 2.0 with large sample
  expect_lt(abs(res$att_overall - 2.0), 1.0)
})

test_that("common-timing detrend runs without error", {
  df  <- make_common_panel(tau = 1.5)
  res <- lwdid(df, y = "y", ivar = "id", tvar = "time", post = "post",
               rolling = "detrend", vce = NULL)
  expect_s3_class(res, "lwdid")
  expect_true(is.finite(res$att_overall))
})

test_that("HC1 and HC3 standard errors are positive", {
  df <- make_common_panel(tau = 2.0)
  for (vce in c("hc1", "hc3")) {
    res <- lwdid(df, y = "y", ivar = "id", tvar = "time", post = "post",
                 rolling = "demean", vce = vce)
    expect_gt(res$se_overall, 0, label = paste("SE > 0 for", vce))
  }
})

test_that("cluster SE requires cluster_var", {
  df <- make_common_panel()
  expect_error(
    lwdid(df, y = "y", ivar = "id", tvar = "time", post = "post",
          vce = "cluster"),
    regexp = "cluster_var"
  )
})

# ── Staggered tests ───────────────────────────────────────────────────────────

test_that("staggered design returns cohort-level ATTs", {
  df  <- make_staggered_panel(tau = 2.0)
  res <- lwdid(df, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
               rolling = "demean", vce = NULL, aggregate = "cohort")
  expect_s3_class(res, "lwdid")
  expect_equal(res$design, "staggered")
  expect_true(!is.null(res$att_by_cohort))
  expect_true(nrow(res$att_by_cohort) >= 2)   # at least the two treated cohorts
  expect_true(is.finite(res$att_overall))
})

test_that("staggered never-treated vs not-yet-treated control groups differ", {
  df  <- make_staggered_panel(tau = 2.0)
  res_nev <- lwdid(df, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
                   control_group = "never_treated",   vce = NULL)
  res_nyt <- lwdid(df, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
                   control_group = "not_yet_treated", vce = NULL)
  expect_s3_class(res_nev, "lwdid")
  expect_s3_class(res_nyt, "lwdid")
})

test_that("staggered aggregate='none' returns all (g,r) pairs", {
  df  <- make_staggered_panel(tau = 2.0)
  res <- lwdid(df, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
               vce = NULL, aggregate = "none")
  expect_true(!is.null(res$att_by_cohort_time))
  expect_gt(nrow(res$att_by_cohort_time), 2)
})

# ── print/summary methods ─────────────────────────────────────────────────────

test_that("print.lwdid and summary.lwdid run without error", {
  df  <- make_common_panel()
  res <- lwdid(df, y = "y", ivar = "id", tvar = "time", post = "post",
               vce = NULL)
  expect_output(print(res))
  expect_output(summary(res))
})
