#' Align two functions
#'
#' This function aligns two functions using SRSF framework. It will align f2
#' to f1
#'
#' @param f1 function 1
#' @param f2 function 2
#' @param time sample points of functions
#' @param lambda controls amount of warping (default = 0)
#' @param pen alignment penalty (default="roughness") options are 
#' second derivative ("roughness"), geodesic distance from id ("geodesic"), and 
#' norm from id ("norm")
#' @param method controls which optimization method (default="DP") options are
#' Dynamic Programming ("DP"), Coordinate Descent ("DP2"), Riemannian BFGS
#' ("RBFGS"), Simultaneous Alignment ("SIMUL"), Dirichlet Bayesian ("dBayes"),
#' and Expo-Map Bayesian ("expBayes")
#' @param w controls LRBFGS (default = 0.01)
#' @param iter number of mcmc iterations for mcmc method (default 2000)
#' @return Returns a list containing \item{f2tilde}{aligned f2}
#' \item{gam}{warping function}
#' @keywords srsf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian
#'  registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @references Lu, Y., Herbei, R., and Kurtek, S. (2017). Bayesian registration
#'  of functions with a Gaussian process prior. Journal of Computational and
#'  Graphical Statistics, DOI: 10.1080/10618600.2017.1336444.
#' @export
#' @examples
#' data("simu_data")
#' out = pair_align_functions(simu_data$f[,1],simu_data$f[,2],simu_data$time)
pair_align_functions <- function(f1, f2, time, lambda=0, pen="roughness",
																 method="DP", w=0.01, iter=2000){

  q1 = f_to_srvf(f1, time)
  q2 = f_to_srvf(f2, time)
  if (method=="dBayes"){
    gam <- pair_align_functions_bayes(f1, f2, time, iter=iter)$gam_a
  } else if (method=="expBayes") {
    gam <- pair_align_functions_expomap(f1, f2, time, iter=iter)$gamma
  } else {
    gam <- optimum.reparam(q1, time, q2, time, lambda, pen, method, w, f1o=f1[1], f2o=f2[1])
  }

  f2_aligned <- warp_f_gamma(f2, time, gam)

  return(list(f2tilde=f2_aligned, gam=gam))
}
