# Event-study plot for the large-N LWDID result. Mirrors the Stata `graph`
# option (WATT(r) by relative period with simultaneous bands; zero/event lines).

#' Event-study plot for a large-N lwdid object
#'
#' Plots WATT(r) against relative time r with simultaneous (sup-t) confidence
#' bands. Requires the \pkg{ggplot2} package.
#'
#' @param x An object of class `"lwdid_large"` from `lwdid(..., method=)`.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot.lwdid_large <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot.lwdid_large requires the 'ggplot2' package.")

  w <- x$watt[!is.na(x$watt$ryear), , drop = FALSE]
  w <- w[order(w$ryear), , drop = FALSE]
  w$phase <- ifelse(w$ryear < 0, "pre", "post")

  ggplot2::ggplot(w, ggplot2::aes(x = ryear, y = watt)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower_ci, ymax = upper_ci,
                   colour = phase), width = 0.15) +
    ggplot2::geom_point(ggplot2::aes(colour = phase), size = 2) +
    ggplot2::geom_line(ggplot2::aes(colour = phase, group = phase)) +
    ggplot2::scale_colour_manual(values = c(pre = "grey40", post = "#2c7fb8"),
                                 guide = "none") +
    ggplot2::labs(
      x = "Time to treatment (r)",
      y = "WATT(r)",
      title = sprintf("LWDID event study (%s, %s)", x$rolling, toupper(x$method)),
      subtitle = "Bars are simultaneous (sup-t) confidence bands") +
    ggplot2::theme_minimal()
}
