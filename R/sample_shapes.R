#' Sample shapes from model
#'
#' @param mu array (n,T) of mean srvf
#' @param K array (2*T,2*T) covariance matrix
#' @param mode Open ("O") or Closed ("C") curves
#' @param no number of principal components
#' @param numSamp number of samples
#' @return Returns a list containing \item{samples}{sample curves}
#' \item{samples.q}{samples srvfs}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'    Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2, parallel=FALSE)
#' # note: use more shapes, small for speed
#' K <- curve_karcher_cov(out$v)
#' out.samples <- sample_shapes(out$mu, K)
sample_shapes <- function(mu, K, mode="O", no=3, numSamp=10){
    n = nrow(mu)
    T1 = ncol(mu)

    out = svd(K)
    U = out$u
    s = out$d
    V = out$v

    if (mode == "O"){
        N = 2
    } else {
        N = 10
    }

    epsilon = 1./(N-1)

    q1 = mu
    q2 = mu
    samples = array(0,dim=c(n,T1,numSamp))
    samples.q = array(0,dim=c(n,T1,numSamp))

    for (i in 1:numSamp){
        v = matrix(0, 2, T1)
        for (m in 1:no){
            v = v + stats::rnorm(1)*sqrt(s[m])*c(U[1:T1,m], U[(T1+1):(2*T1),m])
        }

        q1 = mu
        for (j in 1:(N-1)){
            normv = sqrt(innerprod_q2(v,v))

            if (normv < 1e-4){
                q2 = mu
            } else {
                q2 = cos(epsilon*normv)*q1+sin(epsilon*normv)*v/normv
                if (mode == "C"){
                    q2 = project_curve(q2)
                }
            }

            # Parallel translate tanent vector
            basis2 = find_basis_normal(q2)
            v = parallel_translate(v, q1, q2, basis2, mode)

            q1 = q2
        }

        beta = q_to_curve(q2)
        centroid = calculatecentroid(beta)
        dim(centroid) = c(length(centroid),1)
        beta = beta - repmat(centroid,1,T1)
        samples[,,i] = beta
        samples.q[,,i] = curve_to_q(beta)$q
    }

    return(list(samples=samples,samples.q=samples.q))
}
