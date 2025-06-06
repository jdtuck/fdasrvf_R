#' Align two functions
#'
#' This function aligns the SRVFs of two functions in \eqn{R^1} defined on an
#' interval \eqn{[t_{\min}, t_{\max}]} using dynamic programming or RBFGS
#'
#' @param Q1 A numeric matrix of shape `n_points x n_dimensions` specifying the
#'   SRSF of the 1st `n_dimensions`-dimensional function evaluated on a grid of
#'   size `n_points` of its univariate domain.
#' @param T1 A numeric vector of size `n_points` specifying the grid on which
#'   the 1st SRSF is evaluated.
#' @param Q2 A numeric matrix of shape `n_points x n_dimensions` specifying the
#'   SRSF of the 2nd `n_dimensions`-dimensional function evaluated on a grid of
#'   size `n_points` of its univariate domain.
#' @param T2 A numeric vector of size `n_points` specifying the grid on which
#'   the 1st SRSF is evaluated.
#' @param lambda A numeric value specifying the amount of warping. Defaults to
#'   `0.0`.
#' @param pen alignment penalty (default="roughness") options are
#'   second derivative ("roughness"), geodesic distance from id ("geodesic"), and
#'   norm from id ("l2gam"), srvf norm from id ("l2psi")
#' @param method A string specifying the optimization method. Choices are
#'   `"DP"`, `"DPo"`, `"SIMUL"`, or `"RBFGS"`. Defaults to `"DP"`.
#' @param f1o A numeric vector of size `n_dimensions` specifying the value of
#'   the 1st function at \eqn{t = t_{\min}}. Defaults to `rep(0, n_dimensions)`.
#' @param f2o A numeric vector of size `n_dimensions` specifying the value of
#'   the 2nd function at \eqn{t = t_{\min}}. Defaults to `rep(0, n_dimensions)`.
#' @param nbhd_dim size of the grid (default = 7)
#'
#' @return A numeric vector of size `n_points` storing discrete evaluations of
#'   the estimated boundary-preserving warping diffeomorphism on the initial
#'   grid.
#'
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using Fisher-Rao metric,
#'   arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using phase and amplitude separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#'
#' @export
#' @examples
#' q <- f_to_srvf(simu_data$f, simu_data$time)
#' gam <- optimum.reparam(q[, 1], simu_data$time, q[, 2], simu_data$time)
optimum.reparam <- function(Q1,T1,Q2,T2,
                            lambda = 0,
                            pen = "roughness",
                            method = c("DP", "DPo", "SIMUL", "RBFGS"),
                            f1o = 0.0,
                            f2o = 0.0,
														nbhd_dim=7) {
	pen1 = pen
  pen <- pmatch(pen, c("roughness", "l2gam", "l2psi", "geodesic")) # 1 - roughness, 2 - l2gam, 3 - l2psi, 4 - geodesic
	if (is.na(pen))
    stop("invalid penalty selection")

  M <- length(T1)
  stopifnot(length(T2) == M)

  if (is.null(dim(Q1))) {
    L <- 1
    stopifnot(length(Q1) == M)
    stopifnot(is.null(dim(Q2)) && length(Q2) == M)
  } else {
    L <- nrow(Q1)
    stopifnot(ncol(Q1) == M)
    stopifnot(!is.null(Q2) && nrow(Q2) == L && ncol(Q2) == M)
  }

  method <- match.arg(method, choices = c("DP", "DPo", "SIMUL", "RBFGS"))
  if (method == "DPo" && all(T1 != T2))
    method <- "DP"

  Q1 <- Q1 / pvecnorm(Q1, 2)
  Q2 <- Q2 / pvecnorm(Q2, 2)
  C1 <- srvf_to_f(Q1, T1, f1o)
  C2 <- srvf_to_f(Q2, T2, f2o)

  switch(
    method,
    DP = {
      ret <- DPQ2(Q1, T1, Q2, T2, L, M, M, T1, T2, M, M, lambda,
                  nbhd_dim)
      G <- ret$G[1:ret$size]
      Tf <- ret$T[1:ret$size]
      gam0 <- stats::approx(Tf, G, xout = T2)$y
    },
    DPo = {
      gam0 <- DPQ(Q2, Q1, L, M, lambda, pen, 0)
    },
    SIMUL = {
      if (lambda > 0)
    			warning("penalty not implemented")
      out <- simul_align(C1, C2)
      u <- seq(0, 1, length.out = length(out$g1))
      tmin <- min(T1)
      tmax <- max(T1)
      timet2 <- T1
      timet2 <- (timet2 - tmin) / (tmax - tmin)
      gam0 <- simul_gam(u, out$g1, out$g2, timet2, out$s1, out$s2, timet2)
    },
    RBFGS = {
      time1 <- seq(0, 1, length.out=length(Q1))
      gam0 <- rlbfgs(Q1, Q2, time1, 30, lambda, pen - 1)
      gam0 <- c(gam0)
    }
  )

  (gam0 - gam0[1]) / (gam0[length(gam0)] - gam0[1])  # slight change on scale
}
