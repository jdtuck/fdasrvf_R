#' Place New Functions in Boxplot Regions
#'
#' Aligns new functions to the template of the fit behind an amplitude or
#' phase boxplot (see [align_new()]) and compares their distance to the median
#' with the boxplot regions: the 50% central region, the \eqn{(1-\alpha)}
#' central region, within the whiskers, or outlier.
#'
#' For an `ampbox`, the distance is the amplitude (SRVF) distance of the
#' aligned function to the median SRVF; for a `phbox`, it is the phase distance
#' of its warping function to the median warping function.
#'
#' @param x An object of class `ampbox` or `phbox` typically produced by
#'   [boxplot.fdawarp()] with `what = "stats"`.
#' @param f_new A numeric matrix of shape \eqn{M \times K} of \eqn{K} new
#'   functions sampled at the time grid of the fit, or a numeric vector of
#'   length \eqn{M} for a single function.
#' @param ... Unused here.
#'
#' @return A list with components (one entry or column per function):
#'
#' - `dist`: distance to the median;
#' - `level`: `1` = 50% central region, `2` = \eqn{(1-\alpha)} central region,
#'   `3` = within whiskers, `4` = outlier;
#' - `region`: `level` as text;
#' - `side`: `"Q1"` or `"Q3"`, the quartile direction the function lies in;
#' - `fn`, `qn`, `gam`: the aligned functions, SRVFs and warping functions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' out <- time_warping(simu_data$f, simu_data$time, centroid_type = "median")
#' amp <- boxplot(out, variability_type = "amplitude", what = "stats")
#' classify(amp, simu_data$f[, 1:3])
#' }
classify <- function(x, ...) {
  UseMethod("classify")
}

#' @rdname classify
#' @export
classify.ampbox <- function(x, f_new, ...) {
  check_box_template(x)
  out <- align_to_template(x$template, f_new)
  time <- x$time

  v3 <- x$Q3_q - x$qmedian
  K <- ncol(out$qn)
  dist <- rep(0, K)
  side <- rep(0, K)
  for (k in 1:K) {
    v <- out$qn[, k] - x$qmedian
    dist[k] <- sqrt(trapz(time, v ^ 2))
    side[k] <- trapz(time, v * v3)
  }

  box_classification(x, dist, side, out)
}

#' @rdname classify
#' @export
classify.phbox <- function(x, f_new, ...) {
  check_box_template(x)
  out <- align_to_template(x$template, f_new)

  M <- nrow(out$gam)
  time <- seq(0, 1, length.out = M)
  binsize <- mean(diff(time))
  v3 <- inv_exp_map(x$psi_median, x$Q3_psi)
  K <- ncol(out$gam)
  dist <- rep(0, K)
  side <- rep(0, K)
  for (k in 1:K) {
    psi <- sqrt(pmax(gradient(out$gam[, k], binsize), 0))
    v <- inv_exp_map(x$psi_median, psi)
    dist[k] <- sqrt(trapz(time, v ^ 2))
    side[k] <- trapz(time, v * v3)
  }

  box_classification(x, dist, side, out)
}

check_box_template <- function(x) {
  if (is.null(x$template) || is.null(x$whisker_dis))
    cli::cli_abort(
      "{.arg x} has no region edges; rebuild it with {.code boxplot(..., what = \"stats\")}."
    )
}

box_classification <- function(x, dist, side, aligned) {
  reg <- boxplot_region(dist, c(x$r50, x$ralpha, x$whisker_dis), x$alpha)
  list(
    dist = dist,
    level = reg$level,
    region = reg$region,
    side = ifelse(side >= 0, "Q3", "Q1"),
    fn = aligned$fn,
    qn = aligned$qn,
    gam = aligned$gam
  )
}

# Boxplot region of distances to the median, given the edges
# c(r50, ralpha, whisker_dis): 1 = 50% region, 2 = (1-alpha) region,
# 3 = within whiskers, 4 = outlier
boxplot_region <- function(dist, edges, alpha) {
  names <- c(
    "50% central region",
    sprintf("%g%% central region", 100 * (1 - alpha)),
    "within whiskers",
    "outlier"
  )
  # the whisker can fall inside ralpha, so test it last: the boxplot flags
  # every distance beyond it as an outlier
  level <- rep(3L, length(dist))
  level[dist <= edges[2]] <- 2L
  level[dist <= edges[1]] <- 1L
  level[dist > edges[3]] <- 4L
  list(level = level, region = names[level])
}
