# S3 methods for lwdid class

#' Print an lwdid object
#'
#' Displays a compact summary of estimation results including the design type,
#' transformation method, VCE, overall ATT, SE, t-statistic, and p-value.
#' For staggered designs also prints cohort-specific effects.
#'
#' @param x An object of class `"lwdid"` as returned by `lwdid()`.
#' @param ... Currently ignored.
#' @return Invisibly returns `x`.
#' @export
print.lwdid <- function(x, ...) {
  vce_label <- if (is.null(x$vce)) "OLS" else toupper(x$vce)

  cat("\n")
  cat("Lee-Wooldridge DiD (lwdidR)\n")
  cat(sprintf("Design:      %s\n", x$design))
  cat(sprintf("Transf.:     %s\n", x$rolling))
  cat(sprintf("VCE:         %s\n", vce_label))
  cat(rep("-", 50), "\n", sep = "")

  if (x$design == "staggered") {
    cat(sprintf("Overall ATT: %8.4f\n", x$att_overall))
    cat(sprintf("SE:          %8.4f\n", x$se_overall))
    cat(sprintf("t-stat:      %8.4f\n", x$tstat))
    cat(sprintf("p-value:     %8.4f\n", x$pvalue))
    cat(rep("-", 50), "\n", sep = "")
    cat("\nCohort-specific effects:\n")
    cohort_df <- x$att_by_cohort
    if (!is.null(cohort_df) && nrow(cohort_df) > 0) {
      print(cohort_df[, c("cohort", "att", "se", "tstat", "pvalue")],
            digits = 4, row.names = FALSE)
    }
  } else {
    cat(sprintf("ATT:         %8.4f\n", x$att_overall))
    cat(sprintf("SE:          %8.4f\n", x$se_overall))
    cat(sprintf("t-stat:      %8.4f\n", x$tstat))
    cat(sprintf("p-value:     %8.4f\n", x$pvalue))
    cat(sprintf("95%% CI:      [%.4f, %.4f]\n", x$ci_lower, x$ci_upper))
    cat(sprintf("N (firstpost): %d\n", x$N))
  }
  cat("\n")
  invisible(x)
}


#' Summarise an lwdid object
#'
#' Calls `print.lwdid()` and additionally prints period-specific ATTs
#' (common-timing design) or the first 10 (g, r)-specific ATTs (staggered design).
#'
#' @param object An object of class `"lwdid"` as returned by `lwdid()`.
#' @param ... Currently ignored.
#' @return Invisibly returns `object`.
#' @export
summary.lwdid <- function(object, ...) {
  print(object)
  if (object$design == "staggered" && !is.null(object$att_by_cohort_time)) {
    cat("\n(g, r)-specific effects (first 10 rows):\n")
    gt <- object$att_by_cohort_time
    print(utils::head(gt[, c("cohort", "period", "event_time", "att", "se", "pvalue")], 10),
          digits = 4, row.names = FALSE)
  } else if (!is.null(object$att_by_period)) {
    cat("\nPeriod-specific effects:\n")
    print(object$att_by_period[, c("tindex", "att", "se", "tstat", "pvalue")],
          digits = 4, row.names = FALSE)
  }
  invisible(object)
}
