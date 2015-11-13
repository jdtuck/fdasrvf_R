#' Align two functions
#'
#' This function aligns two srvf functions using Dynamic Programming
#'
#' @param Q1 srvf of function 1
#' @param T1 sample points of function 1
#' @param Q2 srvf of function 2
#' @param T2 sample points of function 2
#' @param lambda controls amount of warping (default = 0)
#' @param method controls which optmization method (default="DP") options are
#' Dynamic Programming ("DP"), Coordinate Descent ("DP2"), Riemannian BFGS
#' ("LRBFGS") and Simultaneous Alignment ("SIMUL")
#' @param w controls LRBFGS (default = 0.01)
#' @param f1o initial value of f1, vector or scalar depending on q1, defaults to zero
#' @param f2o initial value of f2, vector or scalar depending on q1, defaults to zero
#' @return gam warping function
#' @keywords srvf alignment, pca
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(simu_data$f,simu_data$time)
#' gam = optimum.reparam(q[,1],simu_data$time,q[,2],simu_data$time)
optimum.reparam <- function(Q1,T1,Q2,T2,lambda=0,method="DP",w=0.01,f1o=0.0,
                            f2o=0.0){
    n = length(T1)
    Q1=(Q1/pvecnorm(Q1,2))
    Q2=(Q2/pvecnorm(Q2,2))
    G = rep(0,n)
    T = rep(0,n)
    size = 0;
    ret = .Call('DPQ2', PACKAGE = 'fdasrvf', Q1, T1, Q2, T2, 1, n, n, T1, T2, n, n, G, T, size, lambda);

    G = ret$G[1:ret$size]
    Tf = ret$T[1:ret$size]
    gam0 = approx(Tf,G,xout=T2)$y
    gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale
    return(gam)
}
