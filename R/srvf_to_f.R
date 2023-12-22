#' Transformation from SRSF Space
#'
#' This function transforms SRVFs back to the original functional space for
#' functions in \eqn{R^1}.
#'
#' @param q Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the SRSFs that need to be transformed.
#'
#'   - If a vector, it must be of shape \eqn{M} and it is interpreted as a
#'   single \eqn{1}-dimensional curve observed on a grid of size \eqn{M}.
#'   - If a matrix, it must be of shape
#'   \eqn{M \times N}. In this case, it is interpreted as a sample of \eqn{N}
#'   curves observed on a grid of size \eqn{M}, unless \eqn{M = 1} in which case
#'   it is interpreted as a single \eqn{1}-dimensional curve observed on a grid
#'   of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   SRSFs are evaluated.
#' @param f0 Either a numeric value or a numeric vector of or a numeric matrix
#'   specifying the initial value of the curves in the original functional
#'   space. It must be:
#'
#'   - a value if `q` represents a single SRSF.
#'   - a vector of length \eqn{N} if `q` represents a sample of \eqn{N} SRVFs
#'
#' @return A numeric array of the same shape as the input `q` storing the
#'   transformation of the SRVFs `q` back to the original functional space.
#'
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using fisher-rao metric,
#'   arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using amplitude and phase separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#'
#' @export
#' @examples
#' q <- f_to_srvf(simu_data$f, simu_data$time)
#' f <- srvf_to_f(q, simu_data$time, simu_data$f[1, ])
srvf_to_f <- function(q, time, f0 = 0.0) {
  dims <- dim(q)
  dims0 <- dim(f0)
  if (is.null(dims)) { # One uni-dimensional curve
    M <- length(q)
    N <- 1
    integrand <- q * abs(q)
    f <- f0 + cumtrapz(time, integrand)
  } else {
    M <- dims[1]
    N <- dims[2]
    stopifnot(is.null(dims0) && length(f0) == N)
    f <- lapply(1:N, function(n) {
      integrand <- q[, n] * abs(q[, n])
      f0[n] + cumtrapz(time, integrand)
    })
    f <- do.call(cbind, f)
  }
  f
}
