# S3 methods for the large-N LWDID result (class "lwdid_large").

#' Print a large-N lwdid object
#'
#' @param x An object of class `"lwdid_large"` from `lwdid(..., method=)`.
#' @param digits Number of digits for the results table.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.lwdid_large <- function(x, digits = 4, ...) {
  cat("\nLee-Wooldridge DiD (lwdidR) -- large-N\n")
  cat(sprintf("Transf.:  %s    Estimator: %s\n", x$rolling, toupper(x$method)))
  if (!is.null(x$controls))
    cat(sprintf("Controls: %s\n", paste(x$controls, collapse = ", ")))
  cat(sprintf("Controls group: %s    Units: %d    SE: %s\n",
              if (x$never) "never-treated" else "not-yet-treated",
              x$n_units,
              if (isTRUE(x$nose)) "none" else "wild cluster bootstrap"))
  cat(rep("-", 60), "\n", sep = "")

  w <- x$watt
  show_cols <- intersect(c("effect", "ryear", "watt", "se", "t_stat",
                           "p_value", "lower_ci", "upper_ci", "n_cells", "n_units"),
                         names(w))
  print(w[, show_cols], digits = digits, row.names = FALSE)
  cat("\nWATT(r) CIs are simultaneous (sup-t) bands; Pre/Post use pointwise normal CIs.\n")

  if (!is.null(x$attgt)) {
    cat("\nATT(g,t) cells:\n")
    print(x$attgt, digits = digits, row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' Summarise a large-N lwdid object
#' @param object An object of class `"lwdid_large"`.
#' @param ... Ignored.
#' @return Invisibly returns `object`.
#' @export
summary.lwdid_large <- function(object, ...) {
  print(object)
  invisible(object)
}
