#' Interpolate Equally Spaced Points Along a Closed Curve
#'
#' Fits periodic cubic smoothing splines to the supplied coordinates and
#' returns points that are equally spaced by arc length along the fitted curve.
#' The curve is parameterized by cumulative chord length, smoothed using
#' cyclic cubic splines (`mgcv`), and then reparameterized by arc length to
#' obtain equally spaced samples.
#'
#' @param t Either:
#'   \itemize{
#'     \item A single integer specifying the number of equally spaced output
#'       points, or
#'     \item A numeric vector of normalized arc-length positions between 0 and
#'       1 at which to evaluate the curve.
#'   }
#' @param ... Two or more equal-length coordinate vectors (e.g., `x`, `y`, or
#'   `x`, `y`, `z`) describing the input curve.
#' @param n.samples Number of points used to approximate the arc-length
#'   parameterization of the fitted spline. Larger values improve accuracy but
#'   increase computation time.
#' @param k Basis dimension for the cyclic cubic spline. Larger values allow
#'   more flexible fits.
#' @param gamma Smoothing penalty multiplier passed to
#'   `mgcv::gam()`. Values greater than 1 produce smoother curves.
#'
#' @return
#' A matrix with one column per coordinate dimension (`V1`, `V2`, ...),
#' containing interpolated points that are equally spaced by arc length.
#'
#' @details
#' The function first parameterizes the input points by cumulative chord
#' length, fits a cyclic cubic spline to each coordinate using
#' Restricted Maximum Likelihood (REML), estimates arc length from the fitted
#' curve, and then evaluates the splines at parameter values corresponding to
#' equally spaced arc-length positions.
#'
#' @examples
#' x <- cos(seq(0, 2*pi, length.out = 20) * 1.5 - pi/2)
#' y <- sin(seq(0, 2*pi, length.out = 20))
#' plot(x, y, type = "b")
#' lines(interparc(1000, x, y), col = 2)
#' points(interparc(10, x, y), col = 2)
#'
#' @export
interparc <- function(t,
                      ...,
                      n.samples = 5000,
                      k = 20,
                      gamma = 1) {

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required. Install it with:\n",
      "install.packages('mgcv')",
      call. = FALSE
    )
  }

  library(mgcv)

  coords <- list(...)
  n_dim <- length(coords)

  if (n_dim < 2)
    stop("Supply at least x and y coordinates.")

  points <- do.call(cbind, coords)
  n_points <- nrow(points)

  if (any(sapply(coords, length) != n_points))
    stop("All coordinate vectors must have equal length.")

  # If t is an integer, generate equally spaced parameter values.
  if (length(t) == 1 && t > 1 && t %% 1 == 0)
    t <- seq(0, 1, length.out = t)

  ## ----------------------------------------------------
  ## Chord-length parameterization of the original points
  ## ----------------------------------------------------

  segment_lengths <- sqrt(rowSums(diff(points)^2))

  path_param <- c(0, cumsum(segment_lengths))
  path_param <- path_param / max(path_param)

  ## ----------------------------------------------------
  ## Remove duplicated endpoint if the curve is already closed
  ## ----------------------------------------------------

  if (all(points[1, ] == points[n_points, ])) {
    points <- points[-n_points, , drop = FALSE]
    path_param <- path_param[-length(path_param)]
    n_points <- n_points - 1
  }

  ## ----------------------------------------------------
  ## Fit a periodic cubic spline to each coordinate
  ## ----------------------------------------------------

  smooth_fits <- vector("list", n_dim)

  for (j in seq_len(n_dim)) {

    smooth_fits[[j]] <- gam(
      points[, j] ~ s(path_param, bs = "cc", k = k),
      method = "REML",
      gamma = gamma
    )
  }

  ## ----------------------------------------------------
  ## Evaluate the fitted curve to estimate arc length
  ## ----------------------------------------------------

  fit_param <- seq(0, 1, length.out = n.samples)

  fitted_path <- sapply(
    smooth_fits,
    function(fit)
      predict(fit, newdata = data.frame(path_param = fit_param))
  )

  segment_lengths <- sqrt(rowSums(diff(fitted_path)^2))

  arc_length <- c(0, cumsum(segment_lengths))
  arc_length <- arc_length / max(arc_length)

  ## ----------------------------------------------------
  ## Convert desired arc-length positions into spline parameter values
  ## ----------------------------------------------------

  sample_param <- approx(
    x = arc_length,
    y = fit_param,
    xout = t,
    ties = "ordered"
  )$y

  ## ----------------------------------------------------
  ## Evaluate the spline at the equal arc-length locations
  ## ----------------------------------------------------

  result <- sapply(
    smooth_fits,
    function(fit)
      predict(fit, newdata = data.frame(path_param = sample_param))
  )

  colnames(result) <- paste0("V", seq_len(n_dim))

  result
}
