#' Pairwise align two curves
#'
#' This function aligns to curves using Elastic Framework
#'
#' @param beta1 curve 1, provided as a matrix of dimensions \eqn{n \times T} for
#'  \eqn{n}-dimensional curve evaluated on \eqn{T} sample points
#' @param beta2 curve 2, provided as a matrix of dimensions \eqn{n \times T} for
#'  \eqn{n}-dimensional curve evaluated on \eqn{T} sample points
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotation Include rotation (default = `TRUE`)
#' @param scale scale curves to unit length (default = `TRUE`)
#' @return a list containing \item{beta2n}{aligned curve 2 to 1}
#' \item{q2n}{aligned srvf 2 to 1}
#' \item{gam}{warping function}
#' \item{q1}{srvf of curve 1}
#' \item{beta1}{centered curve 1}
#' \item{beta2}{centered curve 2}
#' \item{R}{rotation matrix}
#' \item{tau}{seed}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'    Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_pair_align(beta[, , 1, 1], beta[, , 1, 5])
curve_pair_align <- function(beta1, beta2, mode="O", rotation=TRUE, scale=TRUE){
  T1 = ncol(beta1)
  centroid1 = calculatecentroid(beta1)
  dim(centroid1) = c(length(centroid1),1)
  beta1 = beta1 - repmat(centroid1, 1, T1)
  centroid2 = calculatecentroid(beta2)
  dim(centroid2) = c(length(centroid2),1)
  beta2 = beta2 - repmat(centroid2, 1, T1)

  out = curve_to_q(beta1, scale)
  q1 = out$q
  len1 = out$len
  out = curve_to_q(beta2, scale)
  q2 = out$q
  len2 = out$len
  q2 = curve_to_q(beta2, scale)$q
  out = find_rotation_seed_unique(q1, q2, mode=mode, rotation=rotation, scale=scale)
  gam = out$gambest
  q2n = out$q2best
  R = out$Rbest
  tau = out$tau
  if (mode == "C"){
    beta2n <- shift_f(beta2, tau)
    beta2n[, T1] = beta2n[, 1]
  }
  else
    beta2n = beta2

  beta2n <- R %*% beta2n
  beta2n <- group_action_by_gamma_coord(beta2n, gam)

  ratio <- 1
  if (scale)
    ratio <- len1 / len2

  beta2n = beta2n * ratio

  return(list(beta2n=beta2n, q2n=q2n, gam=gam, q1=q1, beta1=beta1, beta2=beta2,
              R=R, tau=tau))
}
