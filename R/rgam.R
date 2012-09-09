#' Random Warping
#'
#' Generates random warping functions
#'
#' @param N length of warping function
#' @param sigma variance of warping functions
#' @param num number of warping functions
#' @return gam warping functions
#' @keywords diffeomorphism, warping function
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric, 
#'  arXiv:1103.3817v2 [math.ST].
#' @export
#' @examples
#' gam = rgam(N=101, sigma=.01, num=35)
rgam <- function(N, sigma, num){
	library(mvtnorm)
	gam = matrix(0,num,N+1)
	
	time = N + 1
	mu = sqrt(rep(1,N)*time/N)
	
	for (k in 1:num){
		U = runif(N,min=0,max=3)
		C = sigma*diag(U)
		
		v = rmvnorm(1, sigma = C)
		v = v - (mu%*%t(v))%*%mu/time;
		vn = pvecnorm(v,2)/sqrt(time);
		psi = cos(vn)*mu + sin(vn)*v/vn;
		gam[k,] = c(0,cumsum(psi*psi))/time
		
		# numerical stability
		gam[k,] = gam[k,] + (1e-4)*(1:T)/T;
		gam[k,] = (gam[k,] - gam[k,1])/(gam[k,time]-gam[k,1]);
	}           
	
	gam2 = smooth.data(t(gam),10)
	
	gam = t(gam2)
	
	return(gam)
}