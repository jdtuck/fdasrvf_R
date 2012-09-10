#' Convert to SRVF
#'
#' This function converts functions to srvf
#'
#' @param f matrix of functions
#' @param time time
#' @return q matrix of srvfs
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric, 
#'  arXiv:1103.3817v2 [math.ST].
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(f,time)
f_to_srvf <- function(f,time){
	binsize = mean(diff(time))
	eps = .Machine$double.eps
	fy = gradient(f,binsize)
	q = fy/sqrt(abs(fy)+eps)
	return(q)
}