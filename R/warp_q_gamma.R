#' Warp SRSF
#'
#' This function warps srsf \eqn{q} by \eqn{\gamma}
#'
#' @param q vector
#' @param time time
#' @param gamma vector warping function
#' @param spl.int use spline interpolation (default F)
#' @return qnew warped function
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' q <- f_to_srvf(simu_data$f, simu_data$time)
#' qnew <- warp_q_gamma(q[, 1], simu_data$time, seq(0, 1, length.out = 101))
warp_q_gamma <- function(q, time, gamma, spl.int=FALSE){
    M <- length(gamma);
    gam_dev <- gradient(gamma, 1/(M-1))
    if (spl.int){
      qnew <- stats::spline(time,q,xout=(time[length(time)]-time[1])*gamma +
                    time[1])$y*sqrt(gam_dev)
    } else {
      qnew <- stats::approx(time,q,xout=(time[length(time)]-time[1])*gamma +
                    time[1])$y*sqrt(gam_dev)
    }

    return(qnew)
}
