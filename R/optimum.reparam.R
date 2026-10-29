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
#'   no penalty ("none"), second derivative ("roughness"), geodesic distance
#'   from id ("geodesic"), norm from id ("l2gam"), and srvf norm from id
#'   ("l2psi"). The penalty is weighted by `lambda`, so it has no effect when
#'   `lambda = 0`. It is used by the `"DP"`, `"DPo"`, and `"RBFGS"` methods;
#'   `"SIMUL"` ignores it.
#' @param method A string specifying the optimization method. Choices are
#'   `"DP"`, `"DPo"`, `"SIMUL"`, or `"RBFGS"`. Defaults to `"DP"`.
#' @param f1o A numeric vector of size `n_dimensions` specifying the value of
#'   the 1st function at \eqn{t = t_{\min}}. Defaults to `rep(0, n_dimensions)`.
#' @param f2o A numeric vector of size `n_dimensions` specifying the value of
#'   the 2nd function at \eqn{t = t_{\min}}. Defaults to `rep(0, n_dimensions)`.
#' @param nbhd_dim size of the grid for `"DP"` (default = 10, the 63 slopes
#'   `"DPo"` searches)
#' @param scl for `"DPo"` and `"RBFGS"`, the SRVFs are resampled `scl` times
#'   finer to evaluate the alignment cost; the warp stays on the grid `T1`
#'   (default = 3, 1 turns the resampling off). `"DP"` integrates the linearly
#'   interpolated SRVFs exactly, so it does not use `scl`. `"RBFGS"` starts
#'   from a `"DP"` warp on a coarse grid of at most 51 points (using
#'   `nbhd_dim`) and refines it.
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
                            nbhd_dim = 10,
                            scl = 3) {
  # The two solver families number the penalties differently: the DP code
  # (DPQ2, DPQ) takes 0 = none, 1 = roughness, 2 = l2gam, 3 = l2psi,
  # 4 = geodesic, while rlbfgs has no "none" and takes 0 = roughness,
  # 1 = l2gam, 2 = l2psi, 3 = geodesic
  pen <- pmatch(pen, c("none", "roughness", "l2gam", "l2psi", "geodesic"))
	if (is.na(pen))
    stop("invalid penalty selection")
  pen_dp <- pen - 1
  pen_rbfgs <- max(pen_dp - 1, 0)
  # rlbfgs has no "none" code, but a zero weight is the same thing
  lambda_rbfgs <- if (pen_dp == 0) 0 else lambda

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
  if (method == "DPo" && any(T1 != T2))
    method <- "DP"

  Q1 <- Q1 / pvecnorm(Q1, 2)
  Q2 <- Q2 / pvecnorm(Q2, 2)
  C1 <- srvf_to_f(Q1, T1, f1o)
  C2 <- srvf_to_f(Q2, T2, f2o)

  switch(
    method,
    DP = {
      ret <- DPQ2(Q1, T1, Q2, T2, L, M, M, T1, T2, M, M, lambda,
                  nbhd_dim, pen_dp)
      G <- ret$G[1:ret$size]
      Tf <- ret$T[1:ret$size]
      gam0 <- stats::approx(Tf, G, xout = T2)$y
    },
    DPo = {
      gam0 <- DPQ(Q2, Q1, L, M, lambda, pen_dp, 0, scl)
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
      # RBFGS descends from where it starts, so start it from a DP warp on a
      # coarse grid and let it refine that warp on Q1 and Q2 resampled scl
      # times finer, where its linearly interpolated cost is accurate
      time1 <- seq(0, 1, length.out = M)
      nc <- min(M, 51)
      tc <- seq(0, 1, length.out = nc)
      Q1c <- stats::spline(time1, c(Q1), xout = tc)$y
      Q2c <- stats::spline(time1, c(Q2), xout = tc)$y
      ret <- DPQ2(Q1c, tc, Q2c, tc, 1, nc, nc, tc, tc, nc, nc, lambda,
                  nbhd_dim, pen_dp)
      G <- ret$G[1:ret$size]
      Tf <- ret$T[1:ret$size]
      tf <- seq(0, 1, length.out = (M - 1) * scl + 1)
      # start from the exact identity rather than its interpolation
      gamc <- if (identical(G, Tf))
        NULL
      else
        stats::approx(Tf, G, xout = tf)$y
      # lambda is divided by scl to weigh the penalty as it would on time1
      Q1f <- stats::spline(time1, c(Q1), xout = tf)$y
      Q2f <- stats::spline(time1, c(Q2), xout = tf)$y
      lamf <- lambda_rbfgs / scl
      out <- rlbfgs(Q1f, Q2f, tf, 30, lamf, pen_rbfgs, gamc)
      # the DP warp is piecewise linear, which the penalty can make costly to
      # start from, so with a penalty also start from the identity and keep
      # the warp of lower cost
      if (lamf > 0) {
        outi <- rlbfgs(Q1f, Q2f, tf, 30, lamf, pen_rbfgs)
        if (outi$info[3] < out$info[3])
          out <- outi
      }
      # every scl-th point of tf is a point of time1
      gam0 <- c(out$gam)[seq(1, length(tf), by = scl)]
    }
  )

  (gam0 - gam0[1]) / (gam0[length(gam0)] - gam0[1])  # slight change on scale
}

# "norm" is the older name for the "l2gam" penalty, which is the name
# optimum.reparam() expects; any other value is returned unchanged
penalty_alias <- function(pen) {
  if (identical(pen, "norm")) "l2gam" else pen
}
