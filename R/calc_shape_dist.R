#' Elastic Shape Distance
#'
#' Calculates the elastic shape distance between two curves `beta1` and `beta2`.
#' If the input curves are describing multidimensional functional data, then
#' most of the time the user should set `rotation == FALSE`, `scale == FALSE`
#' and `mode == "O"`.
#'
#' @param beta1 A numeric matrix of shape \eqn{L \times M} specifying an
#'   \eqn{L}-dimensional curve evaluated on \eqn{M} sample points.
#' @param beta2 A numeric matrix of shape \eqn{L \times M} specifying an
#'   \eqn{L}-dimensional curve evaluated on \eqn{M} sample points. This curve
#'   will be aligned to `beta1`.
#' @param mode A character string specifying whether the input curves should be
#'   considered open (`mode == "O"`) or closed (`mode == "C"`). Defaults to
#'   `"O"`.
#' @param alignment A boolean value specifying whether the curves should be
#'  aligned before computing the distance matrix. Defaults to `TRUE`.
#' @param rotation A boolean specifying whether the distance should be made
#'   invariant by rotation. Defaults to `TRUE`.
#' @param scale A boolean specifying whether the distance should be made
#'   invariant by scaling. This is effectively achieved by making SRVFs having
#'   unit length and switching to an appropriate metric on the sphere between
#'   normalized SRVFs. Defaults to `TRUE`.
#' @param include.length A boolean specifying whether to include information
#'   about the actual length of the SRVFs in the metric when using normalized
#'   SRVFs to achieve scale invariance. This only applies if `scale == TRUE`.
#'   Defaults to `FALSE`.
#'
#' @return A list with the following components:
#' - `d`: the amplitude (geodesic) distance;
#' - `dx`: the phase distance;
#' - `q1`: the SRVF of `beta1`;
#' - `q2n`: the SRVF of `beta2` after alignment and possible optimal rotation
#' and scaling;
#' - `beta1`: the input curve `beta1`;
#' - `beta2n`: the input curve `beta2` after alignment and possible optimal
#' rotation and scaling.
#' - `R`: the optimal rotation matrix that has been applied to the second curve;
#' - `gam`: the optimal warping function that has been applied to the second
#' curve;
#' - `betascale`: the optimal scaling factor that has been applied to the second
#' curve.
#'
#' @keywords distances
#'
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'   analysis of elastic curves in euclidean spaces. Pattern Analysis and
#'   Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @references Kurtek, S., Srivastava, A., Klassen, E., and Ding, Z. (2012),
#'   “Statistical Modeling of Curves Using Shapes and Related Features,” Journal
#'   of the American Statistical Association, 107, 1152–1165.
#'
#' @export
#' @examples
#' out <- calc_shape_dist(beta[, , 1, 1], beta[, , 1, 4])
calc_shape_dist <- function(beta1, beta2,
                            mode = "O",
                            alignment = TRUE,
                            rotation = TRUE,
                            scale = TRUE,
                            include.length = FALSE) {
  srvf1 <- curve_to_srvf(beta1, scale = scale)
  srvf2 <- curve_to_srvf(beta2, scale = scale)

  out <- match_f2_to_f1(
    srvf1, srvf2, beta2,
    mode = mode,
    alignment = alignment,
    rotation = rotation,
    scale = scale
  )

  list(
    d = out$d,
    dx = out$dx,
    q1 = srvf1$q,
    q2n = out$q2n,
    beta1 = beta1,
    beta2n = out$beta2n,
    R = out$Rbest,
    gam = out$gambest,
    betascale = out$betascale
  )
}
