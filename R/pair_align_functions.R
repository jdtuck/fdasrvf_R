#' Align two functions
#'
#' This function aligns two functions using SRSF framework. It will align f2
#' to f1
#'
#' @param f1 function 1
#' @param f2 function 2
#' @param time sample points of functions
#' @param lambda controls amount of warping (default = 0)
#' @param method controls which optimization method (default="DP") options are
#' Dynamic Programming ("DP"), Coordinate Descent ("DP2"), Riemannian BFGS
#' ("RBFGS") and Simultaneous Alignment ("SIMUL")
#' @param w controls LRBFGS (default = 0.01)
#' @param f1o initial value of f1, vector or scalar depending on q1, defaults to zero
#' @param f2o initial value of f2, vector or scalar depending on q1, defaults to zero
#' @return Returns a list containing \item{f2tilde}{aligned f2}
#' \item{gam}{warping function}
#' @keywords srsf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' out = pair_align_functions(simu_data$f[,1],simu_data$f[,2],simu_data$time)
pair_align_functions <- function(f1, f2, time,lambda=0,method="DP",w=0.01){

  q1 = f_to_srvf(f1, time)
  q2 = f_to_srvf(f2, time)
  gam = optimum.reparam(q1, time, q2, time, lambda, method, w, f1o=f1[1], f2o=f2[1])
  f2_aligned = warp_f_gamma(f2, time, gam)

  return(list(f2tilde=f2_aligned, gam=gam))
}
