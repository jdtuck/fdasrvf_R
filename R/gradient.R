#' Gradient using finite differences
#'
#' This function computes the gradient of `f` using finite differences.
#'
#' @param f Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the curve(s) that need to be differentiated.
#'
#'   - If a vector, it must be of shape \eqn{M} and it is interpreted as a
#'   single \eqn{1}-dimensional curve observed on a grid of size \eqn{M}.
#'   - If a matrix and `multidimensional == FALSE`, it must be of shape
#'   \eqn{M \times N}. In this case, it is interpreted as a sample of \eqn{N}
#'   curves observed on a grid of size \eqn{M}, unless \eqn{M = 1} in which case
#'   it is interpreted as a single \eqn{1}-dimensional curve observed on a grid
#'   of size \eqn{M}.
#'   - If a matrix and `multidimensional == TRUE`,it must be of shape
#'   \eqn{L \times M} and it is interpreted as a single \eqn{L}-dimensional
#'   curve observed on a grid of size \eqn{M}.
#'   - If a 3D array, it must be of shape \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param binsize A numeric value specifying the size of the bins for computing
#'   finite differences.
#' @param multidimensional A boolean specifying if the curves are
#'   multi-dimensional. This is useful when `f` is provided as a matrix to
#'   determine whether it is a single multi-dimensional curve or a collection of
#'   uni-dimensional curves. Defaults to `FALSE`.
#'
#' @return A numeric array of the same shape as the input array `f` storing the
#'   gradient of `f` obtained via finite differences.
#'
#' @keywords srvf alignment
#' @export
#' @examples
#' out <- gradient(simu_data$f[, 1], mean(diff(simu_data$time)))
gradient <- function(f, binsize, multidimensional = FALSE) {
  dims <- dim(f)

  if (is.null(dims)) {
    # f is a single unidimensional curve
    L <- 1
    M <- length(f)
    N <- 1
    g <- rep(0, M)

    # Take forward differences on left and right edges
    g[1] <- (f[2] - f[1]) / binsize
    g[M] <- (f[M] - f[(M - 1)]) / binsize

    # Take centered differences on interior points
    g[2:(M - 1)] <- (f[3:M] - f[1:(M - 2)]) / (2 * binsize)
  } else if (length(dims) == 2) {
    if (multidimensional || dims[1] == 1) {
      # f is single multidimensional curve
      L <- dims[1]
      M <- dims[2]
      N <- 1
      g <- matrix(0, L, M)

      # Take forward differences on left and right edges
      g[, 1] <- (f[, 2] - f[, 1]) / binsize
      g[, M] <- (f[, M] - f[, (M - 1)]) / binsize

      # Take centered differences on interior points
      g[, 2:(M - 1)] <- (f[, 3:M] - f[, 1:(M - 2)]) / (2 * binsize)
    } else {
      # f is multiple unidimensional curves
      L <- 1
      M <- dims[1]
      N <- dims[2]
      g <- matrix(0, M, N)

      # Take forward differences on left and right edges
      g[1, ] <- (f[2, ] - f[1, ]) / binsize
      g[M, ] <- (f[M, ] - f[(M - 1), ]) / binsize

      # Take centered differences on interior points
      g[2:(M - 1), ] <- (f[3:M, ] - f[1:(M - 2), ]) / (2 * binsize)
    }
  }
  else {
    # f is multiple multidimensional curves
    L <- dims[1]
    M <- dims[2]
    N <- dims[3]
    g <- array(0, dim = dims)

    # Take forward differences on left and right edges
    g[, 1, ] <- (f[, 2, ] - f[, 1, ]) / binsize
    g[, M, ] <- (f[, M, ] - f[, (M - 1), ]) / binsize

    # Take centered differences on interior points
    g[, 2:(M - 1), ] <- (f[, 3:M, ] - f[, 1:(M - 2), ]) / (2 * binsize)
  }

  g
}
