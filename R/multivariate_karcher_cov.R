#' Curve Karcher Covariance
#'
#' Calculate Karcher Covariance of a set of curves
#'
#' @param align_data fdacurve object from [multivariate_karcher_mean] of aligned data
#' @return K covariance matrix
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- multivariate_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more shapes, small for speed
#' K <- multivariate_karcher_cov(out)
multivariate_karcher_cov <- function(align_data){
    v = align_data$v
    tmp = dim(v)
    n = tmp[1]
    T = tmp[2]
    N = tmp[3]

    # Compute Karcher covariance of uniformly sampled mean
    N1 = n*T

    tmpv = matrix(0, N1, N)
    for (i in 1:N){
        tmpv[, i] = c(v[, , i])
    }

    K = stats::cov(t(tmpv))

    return(K)
}
