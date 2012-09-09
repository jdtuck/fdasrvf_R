#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding 
#' shooting vectors
#'
#' @param gam matrix (\eqn{M} x \eqn{N}) of \eqn{M} warping functions
#' @return Returns a list containing \item{mu}{mean function}
#' \item{psi}{srvf of warping functions}
#' \item{vec}{shooting vectors}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric, 
#'  arXiv:1103.3817v2 [math.ST].
#' @export
#' @examples
#' out = SqrtMean(gam)
SqrtMean <- function(gam){
	TT = nrow(gam)
	n = ncol(gam)
	
	psi = matrix(0,TT-1,n)
	for (i in 1:n){
		psi[,i] = sqrt(diff(gam[,i])*TT)
	}
	
	# Find Direction
	mu = psi[,1]
	t = 1
	vec = matrix(0,TT-1,n)
	lvm = rep(0,5)
	for (iter in 1:5){
		for (i in 1:n){
			v = psi[,i] - mu
			len = acos(sum(mu*psi[,i])/TT)
			if (len > 0.05){
				vec[,i] = (len/sin(len))*(psi[,i] - cos(len)*mu)
			}else{
				vec[,i] = rep(0,TT-1)
			}
		}
		vm = rowMeans(vec)
		lvm[iter] = sqrt(sum(vm*vm)/TT)
		mu = cos(t*lvm[iter])*mu + (sin(t*lvm[iter])/lvm[iter])*vm
	}
	
	phi = matrix(0,TT,n)
	for (i in 1:n){
		tmp = rep(0,TT)
		tmp[1:TT-1] = psi[,i]*psi[,i]/TT
		phi[,i] = cumsum(tmp)
	}
	
	out = list(mu = mu,psi = psi,vec = vec)
	return(out)
	
}