#' Align two functions
#'
#' This function aligns two srvf functions using Dynamic Programming
#'
#' @param Q1 srvf of function 1
#' @param T1 sample points of function 1
#' @param Q2 srvf of function 2
#' @param T2 sample points of function 2
#' @return gam warping function
#' @keywords srvf alignment, pca
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric, 
#'  arXiv:1103.3817v2 [math.ST].
#' @export
#' @examples
#' gam = optimum.reparam(Q1,T1,Q2,T2)
optimum.reparam <- function(Q1,T1,Q2,T2){
	n = length(time)
	dyn.load(paste("DynamicProgrammingQ2",.Platform$dynlib.ext,sep=""))
	output <-.C("DynamicProgrammingQ2",Q1=as.double(Q1/pvecnorm(Q1,2)),
							T1=as.double(T1),Q2=as.double(Q2/pvecnorm(Q2,2)),T2=as.double(T2),
							m1=as.integer(1),n1=as.integer(n),n2=as.integer(n),tv1=as.double(T1),
							tv2=as.double(T2),n1v=as.integer(n),n2v=as.integer(n),G = as.double(rep(0,n)),
							T = as.double(rep(0,n)),size = as.integer(0))
	G = output$G[1:output$size]
	Tf = output$T[1:output$size]
	gam0 = approx(Tf,G,xout=T2)$y
	gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale
	return(gam)
}