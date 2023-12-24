#' Calculates elastic depth for a set of curves
#'
#' This functions calculates the elastic depth between set of curves. If the
#' curves are describing multidimensional functional data, then
#' `rotated == FALSE` and `mode == 'O'`
#'
#' @param beta Array of sizes \eqn{n \times T \times N} for \eqn{N} curves
#' of dimension \eqn{T} evaluated on a grid of \eqn{n} points
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotated Include rotation (default = `TRUE`)
#' @param scale Include scale (default = `FALSE`)
#' @param parallel run computation in parallel (default = `TRUE`)
#' @return Returns a list containing \item{amp}{amplitude depth}
#' \item{phase}{phase depth}
#' @keywords depth
#' @concept srvf alignment
#' @references T. Harris, J. D. Tucker, B. Li, and L. Shand, "Elastic depths for
#'  detecting shape anomalies in functional data," Technometrics,
#'  10.1080/00401706.2020.1811156, 2020.
#' @export
#' @examples
#' data("mpeg7")
#' # note: use more shapes and iterations, small for speed
#' out = curve_depth(beta[,,1,1:2],maxit=2)
curve_depth <- function(beta, mode="O", rotated=TRUE, scale=FALSE,
                        parallel = FALSE){
  if (parallel){
    cores = max(parallel::detectCores() - 1, 1)
    cl = parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  } else
  {
    foreach::registerDoSEQ()
  }

  dims = dim(beta)
  n <- dims[1]
  T <- dims[2]
  N <- dims[3]

  amp_dist = matrix(0, N, N)
  phs_dist = matrix(0, N, N)
  k = 0

  for (c1 in 1:(N-1)) {

    dist<-foreach::foreach(k = c1:N, .combine=cbind,.packages='fdasrvf') %dopar% {

      out <- calc_shape_dist(beta[,,c1], beta[,,k], mode=mode, rotation=rotated,
                              scale=scale)

      list(out$d,out$dx)
    }

    NN = N-c1+1
    phs = unlist(dist[2,])
    dim(phs)=c(1,NN)
    amp = unlist(dist[1,])
    dim(amp)=c(1,NN)

    phs_dist[c1, c1:N] = phs
    amp_dist[c1, c1:N] = amp
  }

  amp_dist = amp_dist + t(amp_dist)
  phs_dist = phs_dist + t(phs_dist)

  amp = 1 / (1 + apply(amp_dist, 1, stats::median))
  phase = 1 / (1 + apply(phs_dist, 1, stats::median))
  phase = ((2+pi)/pi) * (phase - 2/(2+pi))

  if (parallel) parallel::stopCluster(cl)

  return(list(amp=amp,phase=phase))
}
