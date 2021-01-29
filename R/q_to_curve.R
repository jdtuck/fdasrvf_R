#' Convert to curve space
#'
#' This function converts SRVFs to curves
#'
#' @param q array describing SRVF (n,T)
#' @param scale scale of original beta (default 1)
#' @return beta array describing curve
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' q = curve_to_q(beta[,,1,1])
#' beta1 = q_to_curve(q)
q_to_curve <- function(q, scale=1){
  T1 = ncol(q)
  n = nrow(q)
  qnorm = rep(0, T1)
  for (i in 1:T1) {
    qnorm[i] = pvecnorm(q[, i], 2)
  }
  integrand = matrix(0, n, T1)
  integrand=t(apply(q,1,function(qrow) qrow*qnorm ))
  beta = cumtrapz(1:T1, integrand, 2)/T1

  beta = scale * beta
  
  return(beta)
}
