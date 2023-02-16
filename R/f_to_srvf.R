#' Transformation to SRSF Space
#'
#' This function transforms curves from their original functional space to the
#' SRVF space.
#'
#' @param f Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the functions that need to be transformed.
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
#'   the curves are evaluated.
#' @param multidimensional A boolean specifying if the curves are
#'   multi-dimensional. This is useful when `f` is provided as a matrix to
#'   determine whether it is a single multi-dimensional curve or a collection of
#'   uni-dimensional curves. Defaults to `FALSE`.
#'
#' @return A numeric array of the same shape as the input array `f` storing the
#'   SRSFs of the original curves.
#'
#' @keywords srsf alignment
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
f_to_srvf <- function(f, time, multidimensional = FALSE) {
  binsize <- mean(diff(time))
  eps <- .Machine$double.eps
  g <- gradient(f, binsize, multidimensional = multidimensional)
  # compute norm of g
  dims <- dim(g)
  if (is.null(dims)) {
    # g is a single unidimensional curve (M)
    L <- 1
    M <- length(g)
    N <- 1
    norm_g <- abs(g)
  } else if (length(dims) == 2) {
    if (multidimensional || dims[1] == 1) {
      # g is a single multidimensional curve (LxM)
      L <- dims[1]
      M <- dims[2]
      N <- 1
      norm_g <- matrix(sqrt(colSums(g^2)), nrow = L, ncol = M, byrow = TRUE)
    } else {
      # g is a collection of unidimensional curves (MxN)
      L <- 1
      M <- dims[1]
      N <- dims[2]
      norm_g <- abs(g)
    }
  } else {
    # g is a collection of multidimensional curves (LxMxN)
    L <- dims[1]
    M <- dims[2]
    N <- dims[3]
    norm_g_list <- lapply(1:N, function(n) {
      sqrt(colSums(g[, , n, drop = FALSE]^2))
    })
    norm_g <- array(dim = c(L, M , N))
    for (n in 1:N)
      norm_g[, , n] <- matrix(norm_g_list[[n]], nrow = L, ncol = M, byrow = TRUE)
  }

  g / sqrt(norm_g + eps)
}
