#' Transformation to SRVF Space
#'
#' This function transforms functions in R^1 from their original functional
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
f_to_srvf <- function(f, time) {
  binsize <- mean(diff(time))
  eps <- .Machine$double.eps
  g <- gradient(f, binsize)
  # compute norm of g
  dims <- dim(g)
  if (is.null(dims)) {
    # g is a single unidimensional curve (M)
    L <- 1
    M <- length(g)
    N <- 1
    norm_g <- abs(g)
  } else if (length(dims) == 2) {
      # g is a collection of unidimensional curves (MxN)
      L <- 1
      M <- dims[1]
      N <- dims[2]
      norm_g <- abs(g)
  }

  g / sqrt(norm_g + eps)
}
