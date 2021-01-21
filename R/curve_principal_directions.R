#' Curve PCA
#'
#' Calculate principal directions of a set of curves
#'
#' @param v array (n,T) of shooting vectors
#' @param K array (2*T,2*T) covariance matrix
#' @param mu array (n,T) of mean srvf
#' @param no number of components
#' @param N number of samples on each side of mean
#' @return Returns a list containing \item{s}{singular values}
#' \item{U}{singular vectors}
#' \item{coef}{principal coefficients}
#' \item{pd}{principal directions}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_karcher_mean(beta[,,1,1:2], maxit=2) # note: use more shapes, small for speed
#' K = curve_karcher_cov(out$betamean, out$v)
#' out = curve_principal_directions(out$v, K, out$mu)
curve_principal_directions <- function(v, K, mu, no=3, N=5){
    n = nrow(mu)
    T1 = ncol(mu)

    # SVD
    out = svd(K)
    U = out$u[,1:no]
    s = out$d[1:no]

    # express shapes as coefficients
    tmp = dim(v)
    N = tmp[3]
    VM = apply(v, c(1,2), mean)
    VM = c(VM)
    x = matrix(0, no, N)
    for (ii in 1:N){
        tmpv = v[, , ii]
        x[, ii] = t(U)%*%(c(tmpv)-VM)
    }


    pd = array(list(), c(no, N))
    for (m in 1:no){
        for (i in 1:N){
            tmp = VM + 0.5*(i-5)*sqrt(s[m])*U[,m]
            v1 = tmp
            dim(v1) = c(n,T1)
            q2n = elastic_shooting(mu, v1)
            p = q_to_curve(q2n)

            pd[m, i][[1]] = p

        }
    }

    return(list(s = s, U = U, coef = x, pd = pd, VM=VM))
}
