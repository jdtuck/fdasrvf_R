#' Transformation to SRVF Space
#'
#' This function transforms functions in \eqn{R^1} from their original functional
#' space to the SRVF space.
#'
#' @param f Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the functions that need to be transformed.
#'
#'   - If a vector, it must be of shape \eqn{M} and it is interpreted as a
#'   single \eqn{1}-dimensional curve observed on a grid of size \eqn{M}.
#'   - If a matrix, it must be of shape
#'   \eqn{M \times N}. In this case, it is interpreted as a sample of \eqn{N}
#'   curves observed on a grid of size \eqn{M}, unless \eqn{M = 1} in which case
#'   it is interpreted as a single \eqn{1}-dimensional curve observed on a grid
#'   of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   the functions are evaluated.
#' @param smooth A boolean specifying whether to differentiate a smoothing
#'   spline of `f` instead of the interpolating cubic spline. Defaults to
#'   `FALSE`. With `FALSE`, [srvf_to_f()] inverts this transformation with
#'   \eqn{O(h^4)} accuracy for smooth `f`. With `TRUE`, high-frequency content
#'   is removed on purpose and is not recovered by [srvf_to_f()].
#'
#' @return A numeric array of the same shape as the input array `f` storing the
#'   SRVFs of the original curves.
#'
#' @keywords srvf alignment
#'
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using Fisher-Rao metric,
#'   arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using phase and amplitude Separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#'
#' @export
#' @examples
#' q <- f_to_srvf(simu_data$f, simu_data$time)
f_to_srvf <- function(f, time, smooth = FALSE) {
  eps <- .Machine$double.eps
  if (is.null(dim(f))) {
    g <- spline_derivative(time, f, smooth)
  } else {
    if (length(dim(f)) > 2) {
      stop('wrong input dimensions of f')
    }
    g <- apply(f, 2, function(x) spline_derivative(time, x, smooth))
    if (is.null(dim(g))) g <- matrix(g, nrow = nrow(f))
  }

  g / sqrt(abs(g) + eps)
}

# Derivative of the interpolating cubic spline (or of a smoothing spline when
# `smooth = TRUE`) of `y` evaluated at the knots `time`.
spline_derivative <- function(time, y, smooth = FALSE) {
  if (smooth) {
    sp <- stats::smooth.spline(time, y)
    return(as.numeric(stats::predict(sp, time, deriv = 1)$y))
  }
  as.numeric(stats::splinefun(time, y, method = "fmm")(time, deriv = 1))
}

# Integral from time[1] to each time point of the interpolating cubic spline
# of `y`. Exact inverse of `spline_derivative()` up to the spline error.
spline_cumintegral <- function(time, y) {
  M <- length(time)
  sf <- stats::splinefun(time, y, method = "fmm")
  h <- diff(time)
  idx <- 1:(M - 1)
  b <- sf(time, deriv = 1)[idx]
  c2 <- sf(time, deriv = 2)[idx] / 2
  d <- sf(time, deriv = 3)[idx] / 6
  seg <- h * (y[idx] + h * (b / 2 + h * (c2 / 3 + h * d / 4)))
  c(0, cumsum(seg))
}
