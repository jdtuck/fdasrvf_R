#' Calculates two elastic distance
#'
#' This functions calculates the distances between functions in \eqn{R^1}
#' \eqn{D_y} and \eqn{D_x}, where function 1 is aligned to function 2
#'
#' @param f1 sample function 1, provided as a vector of length \eqn{M}
#' @param f2 sample function 2, provided as a vector of length \eqn{M}
#' @param time sample points of functions, provided as a vector of length
#'    \eqn{M}
#' @param lambda controls amount of warping (default = `0`)
#' @param pen alignment penalty (default = `"roughness"`) options are
#' second derivative (`"roughness"`), geodesic distance from id (`"geodesic"`),
#' and norm from id (`"norm"`)
#' @return Returns a list containing \item{Dy}{amplitude distance}
#' \item{Dx}{phase distance}
#' @keywords distances
#' @concept srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' distances <- elastic.distance(
#'   f1 = simu_data$f[, 1],
#'   f2 = simu_data$f[, 2],
#'   time = simu_data$time
#' )
elastic.distance <- function(f1, f2, time, lambda = 0, pen="roughness"){
    q1 <- f_to_srvf(f1,time)
    q2 <- f_to_srvf(f2,time)
    gam <- optimum.reparam(q1,time,q2,time,lambda,pen)
    fw <- warp_f_gamma(f2,time,gam)
    qw <- warp_q_gamma(q2,time,gam)
    Dy <- sqrt(trapz(time, (q1-qw)^2))

    time1 <- seq(0,1,length.out=length(time))
    binsize <- mean(diff(time1))
    psi <- sqrt(gradient(gam,binsize))
    q1dotq2 = trapz(time1, psi)
    if (q1dotq2 > 1){
      q1dotq2 = 1
    } else if (q1dotq2 < - 1){
      q1dotq2 = -1
    }
    Dx <- acos(q1dotq2)

    return(list(Dy=Dy,Dx=Dx))
}
