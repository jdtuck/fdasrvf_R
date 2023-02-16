#' Transformation from SRSF Space
#'
#' This function transforms SRVFs back to the original functional space.
#'
#' @param q Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the SRSFs that need to be transformed.
#'
#'   - If a vector, it must be of shape \eqn{M} and it is interpreted as a
#'   single \eqn{1}-dimensional curve observed on a grid of size \eqn{M}.
#'   - If a matrix and `multidimensional == FALSE`, it must be of shape
#'   \eqn{M \times N}. In this case, it is interpreted as a sample of \eqn{N}
#'   curves observed on a grid of size \eqn{M}, unless \eqn{M = 1} in which case
#'   it is interpreted as a single \eqn{1}-dimensional curve observed on a grid
#'   of size \eqn{M}.
#'   - If a matrix and `multidimensional == TRUE`, it is interpreted as a single
#'   multi-dimensional curve.
#'   - If a 3D array, it must be of shape \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   SRSFs are evaluated.
#' @param f0 Either a numeric value or a numeric vector of or a numeric matrix
#'   specifying the initial value of the curves in the original functional
#'   space. It must be:
#'
#'   - a value if `q` represents a single \eqn{1}-dimensional SRSF.
#'   - a vector of length \eqn{L} if `q` represents a single
#'   \eqn{L}-dimensional SRSF.
#'   - a vector of length \eqn{N} if `q` represents a sample of \eqn{N}
#'   \eqn{1}-dimensional SRSFs.
#'   - a matrix of shape \eqn{L \times M} if `q` represents a sample of \eqn{N}
#'   \eqn{L}-dimensional  SRSFs.
#' @param multidimensional A boolean specifying if the curves are
#'   multi-dimensional. This is useful when `q` is provided as a matrix to
#'   determine whether it is a single multi-dimensional curve or a collection of
#'   uni-dimensional curves. Defaults to `FALSE`.
#'
#' @return A numeric array of the same shape as the input `q` storing the
#'   transformation of the SRSFs `q` back to the original functional space.
#'
#' @keywords srsf alignment
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
srvf_to_f <- function(q, time, f0 = 0.0, multidimensional = FALSE) {
  dims <- dim(q)
  dims0 <- dim(f0)
  if (is.null(dims)) { # One uni-dimensional curve
    L <- 1
    stopifnot(length(f0) == L)
    M <- length(q)
    N <- 1
    integrand <- q * abs(q)
    f <- f0 + cumtrapz(time, integrand)
  } else {
    if (length(dims) == 2) {
      if (multidimensional || dims[1] == 1) { # One multi-dimensional curve
        L <- dims[1]
        stopifnot(is.null(dims0) && length(f0) == L)
        M <- dims[2]
        N <- 1
        norm_q <- sqrt(colSums(q^2))
        f <- lapply(1:L, function(l) {
          integrand <- q[l, ] * norm_q
          f0[l] + cumtrapz(time, integrand)
        })
        f <- do.call(rbind, f)
      } else { # N uni-dimensional curves
        L <- 1
        M <- dims[1]
        N <- dims[2]
        stopifnot(is.null(dims0) && length(f0) == N)
        f <- lapply(1:N, function(n) {
          integrand <- q[, n] * abs(q[, n])
          f0[n] + cumtrapz(time, integrand)
        })
        f <- do.call(cbind, f)
      }
    } else { # f is 3D array so N multi-dimensional curves and f0 must be (LxN) matrix
      stopifnot(!is.null(dims0) && length(dims0) == 2)
      L <- dims[1]
      stopifnot(dims0[1] == L)
      M <- dims[2]
      N <- dims[3]
      stopifnot(dims0[2] == N)
      norm_q <- lapply(1:N, function(n) {
        sqrt(colSums(q[, , n]^2))
      })
      f <- lapply(1:N, function(n) {
        res <- lapply(1:L, function(l) {
          integrand <- q[l, , n] * norm_q[[n]][l, ]
          f0[l, n] + cumtrapz(time, integrand)
        })
        do.call(rbind, res)
      })
      f <- do.call(cbind, f)
    }
  }
  f
}
