#' Curve PCA
#'
#' Calculate principal directions of a set of curves
#'
#' @param v array of sizes \eqn{n \times T \times N1} for \eqn{N1} shooting
#' vectors of dimension \eqn{n} evaluated on a grid of \eqn{T} points
#' @param K matrix of sizes \eqn{nT \times nT} of covariance matrix
#' @param mu matrix of sizes \eqn{n \times T} of mean srvf
#' @param len length of original curves (default = `NA`)
#' @param no number of components
#' @param N number of samples on each side of mean
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @return Returns a list containing \item{s}{singular values}
#' \item{U}{singular vectors}
#' \item{coef}{principal coefficients}
#' \item{pd}{principal directions}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more shapes, small for speed
#' K <- curve_karcher_cov(out$v)
#' out <- curve_principal_directions(out$v, K, out$mu)
curve_principal_directions <- function(v, K, mu, len = NA, no = 3, N = 5,
                                       mode = "O"){
    n = nrow(mu)
    T1 = ncol(mu)

    # SVD
    out = svd(K)
    U = out$u[,1:no]
    s = out$d[1:no]

    tmp = dim(v)
    N1 = tmp[3]
    VM = apply(v, c(1,2), mean)
    VM = c(VM)
    if (!all(is.na(len))){
        mean_scale = prod(len)^(1/length(len))
        VM = c(VM, mean_scale)
    }

    # express shapes as coefficients
    x = matrix(0, no, N1)
    for (ii in 1:N1){
        tmpv = c(v[, , ii])
        if (!all(is.na(len))){
            tmpv = c(tmpv, len[ii])
        }
        x[, ii] = t(U)%*%(tmpv-VM)
    }

    pd = array(list(), c(no, N))
    for (m in 1:no){
        for (i in 1:N){
            tmp = VM + 0.5*(i-5)*sqrt(s[m])*U[,m]
            if (!all(is.na(len))){
                a = length(tmp)
                v1 = tmp[1:(a-1)]
                tmp_scale = tmp[a]
            } else {
                v1 = tmp
                tmp_scale = 1
            }

            dim(v1) = c(n,T1)
            q2n = elastic_shooting(mu, v1,mode)
            p = q_to_curve(q2n, tmp_scale)

            pd[m, i][[1]] = p

        }
    }

    return(list(s = s, U = U, coef = x, pd = pd, VM=VM))
}
