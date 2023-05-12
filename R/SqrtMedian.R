#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding
#' shooting vectors and finds the median
#'
#' @param gam matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
#' @return Returns a list containing \item{median}{Karcher median psi function}
#' \item{gam_median}{Karcher mean warping function}
#' \item{psi}{srvf of warping functions}
#' \item{vec}{shooting vectors}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' out <- SqrtMedian(simu_warp_median$warping_functions)
SqrtMedian <- function(gam){
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
    vtil <- matrix(0,TT,n)
    d <- rep(0,n)
    dtil <- rep(0,n)
    lvm = rep(0,maxiter)
    iter <- 1

    for (i in 1:n){
      vec[,i] <- inv_exp_map(mu, psi[,i])
      d[i] <- acos(inner_product(mu,psi[,i]))
      vtil[,i] <- vec[,i] / d[i]
      dtil[i] <- 1/d[i]
    }
    vbar <- rowSums(vtil) * sum(dtil)^(-1)
    if (any(is.nan(vbar))){
      vbar <- vec[,1]
    }
    lvm[iter] <- l2_norm(vbar)

    while (lvm[iter]>0.00000001 & iter<maxiter){
      mu <- exp_map(mu, stp*vbar)
      iter <- iter + 1
      for (i in 1:n){
        vec[,i] <- inv_exp_map(mu, psi[,i])
        d[i] <- acos(inner_product(mu,psi[,i]))
        vtil[,i] <- vec[,i] / d[i]
        dtil[i] <- 1/d[i]
      }
      vbar <- rowSums(vtil) * sum(dtil)^(-1)
      lvm[iter] <- l2_norm(vbar)
      if (is.nan(lvm[iter]))
        break
    }

    gam_median = cumtrapz(time, mu*mu)
    gam_median = (gam_median - min(gam_median))/(max(gam_median)-min(gam_median))

    out = list(median = mu,gam_median = gam_median, psi = psi, vec = vec)

    return(out)

}
