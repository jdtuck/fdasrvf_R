#' Curve Karcher Covariance
#'
#' Calculate Karcher Covariance of a set of curves
#'
#' @param v array (n,T,N) for N number of shooting vectors
#' @param len lengths of curves (default=NA)
#' @return K covariance matrix
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more shapes, small for speed
#' K <- curve_karcher_cov(out$v)
curve_karcher_cov <- function(v, len=NA){
    tmp = dim(v)
    M = tmp[1]
    N = tmp[2]
    K = tmp[3]

    # Compute Karcher covariance of uniformly sampled mean
    if (!all(is.na(len))){
        N1 = M*N+1
    } else {
        N1 = M*N
    }
    tmpv = matrix(0,N1,K)
    for (i in 1:K){
        tmp = v[, , i]
        if (!all(is.na(len))){
            tmp = c(tmp, len[i])
        }
        tmpv[, i] = c(tmp)
    }

    K = stats::cov(t(tmpv))

    return(K)
}
