#' Align two curves
#'
#' This function aligns two SRVF functions using Dynamic Programming. If the
#' curves beta1 and beta2 are describing multidimensional functional data, then
#' `rotation == FALSE` and `mode == 'O'`
#'
#' @param beta1 curve 1, provided as a matrix of dimensions \eqn{n \times M} for
#'  \eqn{M}-dimensional curve evaluated on \eqn{n} sample points
#' @param beta2 curve 1, provided as a matrix of dimensions \eqn{n \times M} for
#'  \eqn{M}-dimensional curve evaluated on \eqn{n} sample points
#' @param lambda controls amount of warping (default = `0`)
#' @param method controls which optimization method. Options are
#' Dynamic Programming (`"DP"`). (default = `"DP"`)
#' @param w controls LRBFGS (default = `0.01`)
#' @param rotated boolean if rotation is desired
#' @param isclosed boolean if curve is closed
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @return return a List containing \item{gam}{warping function}
#' \item{R}{rotation matrix}
#' \item{tau}{seed point}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' gam <- reparam_curve(beta[, , 1, 1], beta[, , 1, 5])$gam
reparam_curve <- function(beta1,beta2,lambda=0,method="DP",w=0.01,rotated=T,
                          isclosed=F, mode="O"){
    n1 = nrow(beta2)
    M = ncol(beta2)
    timet = seq(0,1,length.out=M)
    skipm = 4
    auto = 2
    tau = 0
    if (method=="DPo"){
        # Optimize over SO(n) x Gamma
        q1 = curve_to_q(beta1)$q

        # Optimize over SO(n)
        if (rotated){
          out = find_rotation_seed_coord(beta1, beta2, mode)
          beta2 = out$beta2
          R = out$O_hat
          tau = out$tau
        } else{
          R = diag(n1)
          tau = 0
        }
        q2 = curve_to_q(beta2)$q

        # Optimize over Gamma
        q1i = q1
        dim(q1i) = c(M*n1)
        q2i = q2
        dim(q2i) = c(M*n1)
        G = rep(0,M)
        T1 = rep(0,M)
        size = 0
        ret = .Call('DPQ2', PACKAGE = 'fdasrvf', q1i, timet, q2i, timet, n1, M, M, timet, timet, M, M, G, T1, size, lambda, 1);

        G = ret$G[1:ret$size]
        Tf = ret$T[1:ret$size]
        gam0 = stats::approx(Tf,G,xout=timet)$y
    } else if (method=="DP") {
      # Optimize over SO(n) x Gamma
      q1 = curve_to_q(beta1)$q

      # Optimize over SO(n)
      if (rotated){
        out = find_rotation_seed_coord(beta1, beta2);
        beta2 = out$beta2best
        R = out$Rbest
        tau = out$tau
      } else{
        R = diag(n1)
        tau = 0
      }
      q2 = curve_to_q(beta2)$q

      # Optimize over Gamma
      q1 = q1/sqrt(innerprod_q2(q1, q1))
      q2 = q2/sqrt(innerprod_q2(q2, q2))
      q1i = q1
      dim(q1i) = c(M*n1)
      q2i = q2
      dim(q2i) = c(M*n1)
      gam0 = .Call('DPQ', PACKAGE = 'fdasrvf', q1i, q2i, n1, M, lambda, 1, 0, rep(0,M))
    }

    gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale

    return(list(gam=gam,R=R,tau=tau))
}
