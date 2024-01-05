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

  q1 = curve_to_q(beta1, scale)$q
  out = find_rotation_seed_coord(beta1, beta2, mode, rotation, scale)
  gam = out$gambest
  q2n = out$q2best
  beta2n = out$Rbest %*% shift_f(beta2, out$tau)
  beta2n = group_action_by_gamma_coord(beta2n, gam)
  q2n = curve_to_q(beta2n, scale)$q

  return(list(beta2n=beta2n, q2n=q2n, gam=gam, q1=q1, beta1=beta1, beta2=beta2,
              R=out$Rbest, tau=out$tau))
}
