#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding
#' shooting vectors
#'
#' @param gam matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
#' @return Returns a list containing \item{mu}{Karcher mean psi function}
#' \item{gam_mu}{Karcher mean warping function}
#' \item{psi}{srvf of warping functions}
#' \item{vec}{shooting vectors}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_warp")
#' out = SqrtMean(simu_warp$gam)
SqrtMean <- function(gam){
    TT = nrow(gam)
    n = ncol(gam)
    eps = .Machine$double.eps
    time <- seq(0,1,length.out=TT)

    psi = matrix(0,TT,n)
    binsize <- mean(diff(time))
    for (i in 1:n){
        psi[,i] = sqrt(gradient(gam[,i],binsize))
    }

    # Find Direction
    mu = rowMeans(psi)
    stp <- .3
    maxiter = 501
    vec = matrix(0,TT,n)
    lvm = rep(0,maxiter)
    iter <- 1

    for (i in 1:n){
      vec[,i] <- inv_exp_map(mu, psi[,i])
    }
    vbar <- rowMeans(vec)
    lvm[iter] <- l2_norm(vbar)

    while (lvm[iter]>0.00000001 & iter<maxiter){
      mu <- exp_map(mu, stp*vbar)
      iter <- iter + 1
      for (i in 1:n){
        vec[,i] <- inv_exp_map(mu, psi[,i])
      }
      vbar <- rowMeans(vec)
      lvm[iter] <- l2_norm(vbar)
    }

    gam_mu = cumtrapz(time, mu*mu)
    gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))

    out = list(mu = mu,gam_mu = gam_mu, psi = psi, vec = vec)

    return(out)

}
