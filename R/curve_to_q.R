#' Convert to SRVF space
#'
#' This function converts curves or multidimesional functional data to SRVF
#'
#' @param beta a matrix of shape \eqn{n \times T} describing curve or
#'  multidimensional functional data in \eqn{R^n} where \eqn{n} is the dimension
#'  and \eqn{T} is the number of time points
#' @return a numeric array of the same shape as the input array `beta` storing the
#'   SRVFs of the original curves.
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011).
#'  Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and M
#'  achine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' q <- curve_to_q(beta[, , 1, 1])$q
curve_to_q <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    v = apply(beta,1,gradient, 1.0/T1)
    v = t(v)

    q = matrix(0,n,T1)
    len = sqrt(innerprod_q2(v, v))
    for (i in 1:T1){
        L = sqrt(pvecnorm(v[,i],2))
        if (L>0.0001){
            q[,i] = v[,i]/L
        } else {
            q[,i] = v[,i]*0.0001
        }
    }

    len_q = sqrt(innerprod_q2(q, q))
    q = q/len_q

    return(list(q=q,len=len,lenq=len_q))
}
