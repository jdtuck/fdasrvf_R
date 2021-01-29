#' Align Curves
#'
#' Aligns a collection of curves using the elastic square-root velocity (srvf) framework.
#'
#' @param beta array (n,T,N) for N number of curves
#' @param mode Open ("O") or Closed ("C") curves
#' @param rotated Optimize over rotation (default = T)
#' @param maxit maximum number of iterations
#' @return Returns a list containing \item{betan}{aligned curves}
#' \item{qn}{aligned srvfs}
#' \item{betamean}{mean curve}
#' \item{q_mu}{mean SRVFs}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_srvf_align(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
curve_srvf_align <- function(beta, mode="O", rotated=T, maxit=20,ms = "mean"){
    if (mode=="C"){
      isclosed = TRUE
    }
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]
    for (i in 1:N){
      beta[,,i] = beta[,,i]/sqrt(innerprod_q2(beta[,,i],beta[,,i]))
      beta1 = beta[,,i]
      centroid1 = calculatecentroid(beta1)
      dim(centroid1) = c(length(centroid1),1)
      beta[,,i] = beta1 - repmat(centroid1,1,T1)
    }
    out = curve_karcher_mean(beta, mode, rotated, maxit,ms)
    mu = out$mu
    betamean = out$betamean/sqrt(innerprod_q2(out$betamean,out$betamean))
    v = out$v
    q = out$q

    qn = array(0, c(n,T1,N))
    betan = array(0, c(n,T1,N))
    
    # align to mean
    for (ii in 1:N){
        q1 = q[,,ii]
        beta1 = beta[,,ii]

        out = find_rotation_seed_unqiue(mu,q1,mode)
        beta1 = out$Rbest%*%beta1
        beta1n = group_action_by_gamma_coord(beta1, out$gambest)
        q1n = curve_to_q(beta1n)

        out = find_best_rotation(mu, q1n)
        qn[,,ii] = out$q2new
        btmp = out$R%*%beta1n
        betan[,,ii] = btmp/sqrt(innerprod_q2(btmp,btmp))

    }

    return(list(betan=betan, qn=qn, betamean=betamean, q_mu=mu))
}
