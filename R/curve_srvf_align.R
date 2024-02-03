#' Align Curves
#'
#' Aligns a collection of curves using the elastic square-root velocity (srvf)
#' framework. If the curves are describing multidimensional functional data, then
#' `rotated == FALSE` and `mode == 'O'`
#'
#' @param beta Array of sizes \eqn{n \times T \times N} for \eqn{N} curves
#' of dimension \eqn{n} evaluated on a grid of \eqn{T} points
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotated Optimize over rotation (default = `TRUE`)
#' @param scale scale curves to unit length (default = `TRUE`)
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit maximum number of iterations
#' @param ms string defining whether the Karcher mean ("mean") or Karcher median
#'   ("median") is returned (default = "mean")
#' @param parallel A boolean specifying whether to run calculations in parallel.
#'   Defaults to `TRUE`.
#' @return Returns a list containing \item{betan}{aligned curves}
#' \item{qn}{aligned srvfs}
#' \item{betamean}{mean curve}
#' \item{q_mu}{mean SRVFs}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and
#'    Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' # note: use more shapes and iterations, small for speed
#' out = curve_srvf_align(beta[,,1,1:2],maxit=2,parallel=FALSE)
curve_srvf_align <- function(beta, mode = "O", rotated = TRUE, scale = TRUE,
                             lambda = 0.0, maxit = 20, ms = "mean",
                             parallel=TRUE){

    if (mode == "C"){
      isclosed = TRUE
    }
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]

    for (ii in 1:N) {
      beta1 = beta[ , , ii]
      centroid1 = calculatecentroid(beta1)
      dim(centroid1) = c(length(centroid1), 1)
      beta1 = beta1 - repmat(centroid1, 1, T1)
      beta[ , , ii] = beta1
    }

    out = curve_karcher_mean(beta, mode, rotated, scale, lambda, maxit, ms, parallel)
    beta<-out$beta
    mu = out$mu
    betamean = out$betamean
    v = out$v
    q = out$q

    qn = array(0, c(n, T1, N))
    betan = array(0, c(n, T1, N))
    rotmat = array(0, c(n, n, N))
    gams = matrix(0, T1, N)

    if (parallel) {
      cores <- max(parallel::detectCores() - 1, 1)
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
    } else
      foreach::registerDoSEQ()

    # align to mean
    outfor <- foreach::foreach(n = 1:N, .combine = cbind, .packages='fdasrvf') %dopar% {
      out <- curve_align_sub(beta[, , n], q[, , n], mu, mode, rotated, scale, lambda)

      list(out$qn, out$betan, out$rotmat, out$gam)
    }

    qn <- unlist(outfor[1, ])
    dim(qn) <- c(n, T1, N)

    betan <- unlist(outfor[2, ])
    dim(betan) <- c(n, T1, N)

    rotmat <- unlist(outfor[3, ])
    dim(rotmat) <- c(n, n, N)

    gams <- unlist(outfor[4, ])
    dim(gams) <- c(T1, N)

    if (parallel) parallel::stopCluster(cl)

    return(list(betan = betan, qn = qn, betamean = betamean, q_mu = mu,
                rotmat = rotmat, gams = gams, v = v))
}
