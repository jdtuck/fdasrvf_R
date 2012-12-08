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
#'  @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Dec 2012. Generative Models for Function Data using Phase and Amplitude Separation, 
#'  accepted to Computational Statistics and Data Analysis.
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(simu_data$f,simu_data$time)
f_to_srvf <- function(f,time){
	binsize = mean(diff(time))
	eps = .Machine$double.eps
	tmp = gradient.spline(f,binsize)
	q = tmp$g/sqrt(abs(tmp$g)+eps)
	return(q)
}