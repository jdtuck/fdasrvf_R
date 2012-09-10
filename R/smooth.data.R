#' Smooth Functions
#'
#' This function smooths functions using standard box filter
#'
#' @param f matrix (\eqn{M} x \eqn{N}) of \eqn{M} functions with \eqn{N} samples 
#' @param sparam number of times to run box filter
#' @return fo smoothed functions
#' @keywords srvf alignment, pca
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric, 
#'  arXiv:1103.3817v2 [math.ST].
#' @export
#' @examples
#' data("simu_data")
#' fo = smooth.data(simu_data$f,25)
smooth.data <- function(f,sparam){
	M = nrow(f)
	N = ncol(f)
	for (r in 1:sparam){
		for (i in 1:N){
			f[2:(M-1),i] = (f[1:(M-2),i]+2*f[2:(M-1),i] + f[3:M,i])/4
		}
	}
	return(f)
}