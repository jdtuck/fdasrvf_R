#' Warp Function
#'
#' This function warps function \eqn{f} by \eqn{\gamma}
#'
#' @param f vector function
#' @param time time
#' @param gamma vector warping function
#' @param spl.int use spline interpolation (default F)
#' @return fnew warped function
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' fnew <- warp_f_gamma(
#'   f = simu_data$f[, 1],
#'   time = simu_data$time,
#'   gamma = seq(0, 1, length.out = 101)
#' )
warp_f_gamma <- function(f, time, gamma, spl.int=FALSE){
    if (spl.int){
      fnew <- spline(time,f,xout=(time[length(time)]-time[1])*gamma + time[1])$y
    } else {
      fnew <- approx(time,f,xout=(time[length(time)]-time[1])*gamma + time[1])$y
    }

    return(fnew)
}
