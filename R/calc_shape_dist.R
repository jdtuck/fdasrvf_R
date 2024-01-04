#' Elastic Shape Distance
#'
#' Calculate elastic shape distance between two curves beta1 and beta2. If the
#' curves beta1 and beta2 are describing multidimensional functional data, then
#' `rotation == FALSE` and `mode == 'O'`
#'
#' @param beta1 curve1, provided as a matrix of sizes \eqn{n\times T} for
#'  \eqn{n}-dimensional curve on \eqn{T} sample points
#' @param beta2 curve 2, provided as a matrix of sizes \eqn{n\times T} for
#'  \eqn{n}-dimensional curve on \eqn{T} sample points
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotation Include rotation (default = `TRUE`)
#' @param scale scale curves to unit length (default = `TRUE`)
#' @param include.length include length in distance calculation (default = `FALSE`)
#'                       this only applies if `scale=TRUE`
#' @return Returns a list containing \item{d}{geodesic distance}
#' \item{dx}{phase distance}
#' \item{q1}{srvf of curve 1}
#' \item{q2n}{srvf of aligned curve 2}
#' @keywords distances
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'  analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'  Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @references Kurtek, S., Srivastava, A., Klassen, E., and Ding, Z. (2012),
#'  “Statistical Modeling of Curves Using Shapes and Related Features,” Journal
#'  of the American Statistical Association, 107, 1152–1165.
#' @export
#' @examples
#' out <- calc_shape_dist(beta[, , 1, 1], beta[, , 1, 4])
calc_shape_dist <- function(beta1, beta2, mode="O", rotation=TRUE,
                            scale=TRUE, include.length=FALSE){
  T1 = ncol(beta1)
  centroid1 = calculatecentroid(beta1)
  dim(centroid1) = c(length(centroid1),1)
  beta1 = beta1 - repmat(centroid1,1,T1)
  out1 = curve_to_q(beta1, scale)
  q1 = out1$q
  lenq1 = out1$lenq
  lenq2 = curve_to_q(beta2, scale)$lenq

  centroid1 = calculatecentroid(beta2)
  dim(centroid1) = c(length(centroid1),1)
  beta2 = beta2 - repmat(centroid1,1,T1)

  out = find_rotation_seed_coord(beta1, beta2, mode, rotation, scale)
  q1dotq2 = innerprod_q2(q1, out$q2best)

  if (q1dotq2 > 1){
    q1dotq2 = 1
  } else if(q1dotq2 < -1){
    q1dotq2 = -1
  }

  if (scale){
    if (include.length){
      d = sqrt(acos(q1dotq2)^2+log(lenq1/lenq2)^2)
    } else{
      d = acos(q1dotq2)
    }
  } else {
    d = sqrt(sum((q1-out$q2best)^2)/T1)
  }

  gam = out$gambest
  time1 <- seq(0,1,length.out=T1)
  binsize <- mean(diff(time1))
  psi <- sqrt(gradient(gam,binsize))
  v <- inv_exp_map(rep(1,length(gam)), psi)
  dx <- sqrt(trapz(time1, v^2))

  return(list(d=d,dx=dx,q1=q1,q2n=out$q2best))
}
