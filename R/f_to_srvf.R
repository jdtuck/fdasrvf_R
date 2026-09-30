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
    dims <- dim(f)
    if (length(dims) > 2 && dims[1] > 1) {
      stop('wrong input dimensions of f')
    }
    # a 1 x M x N array is treated as an M x N matrix
    fm <- if (length(dims) > 2) matrix(f, dims[2], dims[3]) else f
    g <- apply(fm, 2, function(x) spline_derivative(time, x, smooth))
    g <- array(g, dim = dims)
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

# Integral from time[1] to each time point of `y`, using the trapezoid rule
# corrected with the slopes of the interpolating cubic spline of `y`
# (Hermite/Euler-Maclaurin end correction). This is O(h^4) accurate and
# inverts `spline_derivative()` more tightly than either the plain trapezoid
# rule or integrating the spline itself.
spline_cumintegral <- function(time, y) {
  M <- length(time)
  h <- diff(time)
  s <- stats::splinefun(time, y, method = "fmm")(time, deriv = 1)
  c(0, cumsum(h * (y[-1] + y[-M]) / 2 + h^2 * (s[-M] - s[-1]) / 12))
}
