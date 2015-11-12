#' Align two functions
#'
#' This function aligns two srvf functions using Dynamic Programming
#'
#' @param Q1 srvf of function 1
#' @param T1 sample points of function 1
#' @param Q2 srvf of function 2
#' @param T2 sample points of function 2
#' @param lambda controls amount of warping (default = 0)
#' @return gam warping function
#' @keywords srvf alignment, pca
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#'  @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(simu_data$f,simu_data$time)
#' gam = optimum.reparam(q[,1],simu_data$time,q[,2],simu_data$time)
optimum.reparam <- function(Q1,T1,Q2,T2,lambda = 0){
    n = length(T1)
    output <-.C("DynamicProgrammingQ2",Q1=as.double(Q1/pvecnorm(Q1,2)),
                            T1=as.double(T1),Q2=as.double(Q2/pvecnorm(Q2,2)),T2=as.double(T2),
                            m1=as.integer(1),n1=as.integer(n),n2=as.integer(n),tv1=as.double(T1),
                            tv2=as.double(T2),n1v=as.integer(n),n2v=as.integer(n),G = as.double(rep(0,n)),
                            T = as.double(rep(0,n)),size = as.integer(0),lam1 = as.double(lambda))
    G = output$G[1:output$size]
    Tf = output$T[1:output$size]
    gam0 = approx(Tf,G,xout=T2)$y
    gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale
    return(gam)
}
