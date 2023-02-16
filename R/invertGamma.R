#' Invert Warping Function
#'
#' This function calculates the inverse of gamma
#'
#' @param gam vector of \eqn{N} samples
#' @return Returns gamI inverted vector
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' out <- invertGamma(simu_warp$gam[, 1])
invertGamma <- function(gam){
    N = length(gam)
    x = (0:(N-1))/(N-1)
    gamI = approx(gam,x,xout=x)$y
    gamI[N] = 1
    gamI = gamI/gamI[N]
    return(gamI)
}
