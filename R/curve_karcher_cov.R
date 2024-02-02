#' Curve Karcher Covariance
#'
#' Calculate Karcher Covariance of a set of curves
#'
#' @param v array of sizes \eqn{n \times T \times N} for \eqn{N} shooting
#' vectors of dimension \eqn{n} evaluated on a grid of \eqn{T} points
#' @param len lengths of curves (default = `NA`)
#' @return K covariance matrix
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2, parallel=FALSE)
#' # note: use more shapes, small for speed
#' K <- curve_karcher_cov(out$v)
curve_karcher_cov <- function(v, len = NA){
    tmp = dim(v)
    n = tmp[1]
    T = tmp[2]
    N = tmp[3]

    # Compute Karcher covariance of uniformly sampled mean
    if (!all(is.na(len))){
        N1 = n*T + 1
    } else {
        N1 = n*T
    }
    tmpv = matrix(0, N1, N)
    for (i in 1:N){
        tmp = v[, , i]
        if (!all(is.na(len))){
            tmp = c(tmp, len[i])
        }
        tmpv[, i] = c(tmp)
    }

    K = stats::cov(t(tmpv))

    return(K)
}
