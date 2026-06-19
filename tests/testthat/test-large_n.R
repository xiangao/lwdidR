# Tests for the large-N path (RA / IPW / IPWRA, WATT(r), inference).

# ---- helper: synthetic staggered panel with a known dynamic effect ----
make_panel <- function(N = 150, seed = 1) {
  set.seed(seed)
  yrs <- 2000:2010
  g_opts <- c(0, 2004, 2006, 2008)
  gid <- sample(g_opts, N, replace = TRUE, prob = c(.4, .2, .2, .2))
  do.call(rbind, lapply(seq_len(N), function(i) {
    ui <- rnorm(1); x1 <- rnorm(1); x2 <- runif(1)
    data.frame(id = i, year = yrs, g = gid[i], x1 = x1, x2 = x2,
      y = ui + 0.05 * (yrs - 2000) + rnorm(length(yrs), 0, 0.3) +
          ifelse(gid[i] > 0 & yrs >= gid[i], 0.4 * (yrs - gid[i] + 1) + 0.1 * x1, 0))
  }))
}

test_that("large-N recovers a known dynamic effect (all estimators)", {
  df <- make_panel()
  for (m in c("ra", "ipw", "ipwra")) {
    res <- lwdid(df, "y", "id", "year", gvar = "g", rolling = "demean",
                 method = m, controls = c("x1", "x2"), reps = 199, seed = 1)
    expect_s3_class(res, "lwdid_large")
    w <- res$watt
    # WATT(r) should track 0.4*(r+1) for r = 0,1,2
    for (rr in 0:2) {
      got <- w$watt[w$effect == sprintf("WATT(%d)", rr)]
      expect_equal(got, 0.4 * (rr + 1), tolerance = 0.1)
    }
    # demean anchors r = -1 to exactly 0
    expect_equal(w$watt[w$effect == "WATT(-1)"], 0)
    expect_equal(w$se[w$effect == "WATT(-1)"], 0)
    # bootstrap SE should be close to the analytic cluster-robust SE
    post <- w[w$effect == "Post_avg", ]
    expect_equal(post$se, post$se_analytic, tolerance = 0.15)
  }
})

test_that("detrend anchors r in {-1,-2}", {
  df <- make_panel()
  res <- lwdid(df, "y", "id", "year", gvar = "g", rolling = "detrend",
               method = "ra", reps = 99, seed = 1)
  w <- res$watt
  expect_equal(w$watt[w$effect == "WATT(-1)"], 0)
  expect_equal(w$watt[w$effect == "WATT(-2)"], 0)
})

test_that("ipw/ipwra require controls", {
  df <- make_panel()
  expect_error(lwdid(df, "y", "id", "year", gvar = "g", method = "ipw"),
               "controls")
})

# ---- oracle validation against Stata lwdid v2.4.2 (skips until do-file is run) ----
ref_dir <- normalizePath(
  file.path(testthat::test_path(), "..", "..", "data-raw", "stata-reference"),
  mustWork = FALSE)

# map a config label to lwdid() arguments
cfg_args <- function(label) {
  a <- list(rolling = if (grepl("detrend", label)) "detrend" else "demean",
            method  = if (grepl("ipwra", label)) "ipwra" else if (grepl("ipw", label)) "ipw" else "ra",
            controls = c("x1", "x2", "x3"), pre = -1L, never = FALSE)
  if (grepl("_nox",  label)) a$controls <- NULL
  if (grepl("_pre3", label)) a$pre <- 3L
  if (grepl("_never", label)) a$never <- TRUE
  a
}

test_that("large-N WATT(r) matches Stata lwdid on lw_walmart", {
  walmart_csv <- file.path(ref_dir, "lw_walmart.csv")
  refs <- Sys.glob(file.path(ref_dir, "ref_walmart_*.csv"))
  skip_if(!file.exists(walmart_csv) || length(refs) == 0,
          "Stata reference CSVs not generated yet (run data-raw/stata-reference/lwdid_reference.do).")

  d <- utils::read.csv(walmart_csv)
  for (rf in refs) {
    label <- sub("^ref_walmart_(.*)\\.csv$", "\\1", basename(rf))
    a <- cfg_args(label)
    res <- lwdid(d, "log_wholesale_emp", "cid", "year", gvar = "first_year",
                 rolling = a$rolling, method = a$method, controls = a$controls,
                 pre = a$pre, never = a$never, reps = 99, seed = 1)
    ref <- utils::read.csv(rf)
    rr <- res$watt
    # WATT(r) point estimates: deterministic, must match Stata to ~machine eps
    wr <- merge(rr[!is.na(rr$ryear), c("ryear", "watt")],
                ref[!is.na(ref$ryear), c("ryear", "watt")],
                by = "ryear", suffixes = c("_r", "_stata"))
    expect_equal(wr$watt_r, wr$watt_stata, tolerance = 1e-5,
                 info = paste("WATT(r) point estimates,", label))

    # Standard errors: compare the deterministic analytic cluster-robust SE to
    # Stata's wild-bootstrap SE (by effect label, incl. Pre/Post). They agree up
    # to bootstrap Monte-Carlo error, so check the median relative gap.
    sem <- merge(rr[, c("effect", "se_analytic")], ref[, c("effect", "se")],
                 by = "effect")
    sem <- sem[is.finite(sem$se) & sem$se > 0, ]
    rel <- abs(sem$se_analytic - sem$se) / sem$se
    expect_lt(stats::median(rel), 0.06,
              label = paste("median rel. SE gap,", label))
  }
})

test_that("small-N path matches Stata lwdid on lw_smoking", {
  smoking_csv <- file.path(ref_dir, "lw_smoking.csv")
  ref_csv <- file.path(ref_dir, "ref_smoking_smallN.csv")
  skip_if(!file.exists(smoking_csv) || !file.exists(ref_csv),
          "Stata reference CSVs not generated yet.")
  d <- utils::read.csv(smoking_csv)
  ref <- utils::read.csv(ref_csv)
  d$gvarNA <- ifelse(d$first_year == 0, NA, d$first_year)
  for (roll in c("demean", "detrend")) {
    res <- lwdid(d, "lcigsale", "state", "year", gvar = "gvarNA", rolling = roll)
    att_stata <- ref$value[ref$config == roll & ref$key == "att"]
    expect_equal(res$att_overall, att_stata, tolerance = 1e-4,
                 info = paste("small-N ATT,", roll))
  }
})
