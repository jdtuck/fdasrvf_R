# Helpers shared by the Bayesian infinite-dimensional HMC alignment routines
# (`pair_align_functions_infHMC()` and `time_warping_bayes_infHMC()`).
#
# Ported from the MATLAB package (util/bayes_warp_basis.m,
# util/bayes_coef_to_gam.m, util/bayes_gam_to_coef.m, util/bayes_warp_geom.m,
# util/infhmc_step.m, util/infhmc_step_batch.m).

#' Zero-mean orthonormal basis for warping coefficients
#'
#' Returns `ncoef` basis functions on `t` in \eqn{[0,1]} that integrate to
#' zero (so that they lie in the tangent space at the identity / give a
#' centered \eqn{h}) and are orthonormal in \eqn{L^2[0,1]}.
#'
#' @param t grid on \eqn{[0,1]}
#' @param ncoef number of basis functions (must be even for `"fourier"`)
#' @param type either `"fourier"` (\eqn{\sqrt 2 \sin(2 \pi j t)},
#'   \eqn{\sqrt 2 \cos(2 \pi j t)}) or `"legendre"` (shifted Legendre
#'   polynomials of degree 1 to `ncoef`)
#'
#' @return A `length(t)` by `ncoef` basis matrix.
#'
#' @keywords internal
bayes_warp_basis <- function(t, ncoef, type = c("fourier", "legendre")) {
  type <- match.arg(tolower(type[1]), c("fourier", "legendre"))
  t <- as.numeric(t)
  Phi <- matrix(0, length(t), ncoef)
  if (type == "fourier") {
    for (i in 1:ncoef) {
      j <- ceiling(i / 2)
      if (i %% 2 == 1)
        Phi[, i] <- sqrt(2) * sin(2 * j * pi * t)
      else
        Phi[, i] <- sqrt(2) * cos(2 * j * pi * t)
    }
  } else {
    # shifted Legendre polynomials P_i(2t - 1), normalized in L2[0,1]
    x <- 2 * t - 1
    pm1 <- rep(1, length(x))
    p <- x
    for (i in 1:ncoef) {
      if (i > 1) {
        pnew <- ((2 * i - 1) * x * p - (i - 1) * pm1) / i
        pm1 <- p
        p <- pnew
      }
      Phi[, i] <- sqrt(2 * i + 1) * p
    }
  }
  Phi
}

#' Warping functions from basis coefficients
#'
#' @param c a \eqn{D} by \eqn{K} matrix of coefficients, one column per
#'   warping function
#' @param Phi a \eqn{T} by \eqn{D} basis matrix on the uniform grid
#'   `seq(0, 1, length.out = T)`
#' @param rep `0` for the shooting vector (v-transform), `1` for the clr
#'   (h-transform) representation
#'
#' @return A list with the \eqn{T} by \eqn{K} matrices `gam` of warping
#'   functions and `psi` of their SRSFs \eqn{\sqrt{\gamma'}}.
#'
#' @keywords internal
bayes_coef_to_gam <- function(c, Phi, rep) {
  c <- as.matrix(c)
  TT <- nrow(Phi)
  dt <- 1 / (TT - 1)
  V <- Phi %*% c
  K <- ncol(c)
  psi <- matrix(0, TT, K)
  for (k in 1:K) {
    v <- V[, k]
    if (rep == 0) {
      r <- sqrt(sum(v ^ 2) * dt)
      if (r < 1e-10)
        p <- rep(1, TT) + v
      else
        p <- cos(r) + sin(r) / r * v
      psi[, k] <- abs(p)
    } else {
      psi[, k] <- exp((v - max(v)) / 2)
    }
  }
  psi2 <- psi ^ 2
  G <- rbind(rep(0, K), apply(dt / 2 * (psi2[-TT, , drop = FALSE] +
                                          psi2[-1, , drop = FALSE]), 2, cumsum))
  GT <- G[TT, ]
  gam <- sweep(G, 2, GT, "/")
  gam[TT, ] <- 1
  psi <- sweep(psi, 2, sqrt(GT), "/")
  list(gam = gam, psi = psi)
}

#' Project warping functions onto warping basis coefficients
#'
#' Inverse of [bayes_coef_to_gam()] up to basis truncation.
#'
#' @param gam a \eqn{T} by \eqn{K} matrix of warping functions on the uniform
#'   grid `seq(0, 1, length.out = T)`
#' @param Phi an orthonormal zero-mean basis, see [bayes_warp_basis()]
#' @param rep `0` for the shooting vector (v-transform), `1` for the clr
#'   (h-transform) representation
#'
#' @return A \eqn{D} by \eqn{K} matrix of coefficients.
#'
#' @keywords internal
bayes_gam_to_coef <- function(gam, Phi, rep) {
  gam <- as.matrix(gam)
  TT <- nrow(gam)
  K <- ncol(gam)
  dt <- 1 / (TT - 1)
  w <- rep(dt, TT)
  w[c(1, TT)] <- dt / 2
  cc <- matrix(0, ncol(Phi), K)
  for (k in 1:K) {
    g <- (gam[, k] - gam[1, k]) / (gam[TT, k] - gam[1, k])
    gd <- pmax(gradient(g, dt), 1e-8)
    if (rep == 0) {
      psi <- sqrt(gd)
      psi <- psi / sqrt(sum(w * psi ^ 2))
      v <- inv_exp_map(rep(1, TT), psi)
      r <- sqrt(sum(w * v ^ 2))
      if (r >= 0.95 * pi)
        v <- v * (0.95 * pi / r)
    } else {
      v <- log(gd)
      v <- v - sum(w * v)
    }
    cc[, k] <- t(Phi) %*% (w * v)
  }
  cc
}

#' Negative log-likelihood and exact gradient of a warp
#'
#' Model in SRVF space: \eqn{q_1(t_k) = (q_2 \circ \gamma)(t_k)
#' \sqrt{\gamma'(t_k)} + e_k}, \eqn{e_k \sim N(0, \mathrm{var}_1)}, with
#' \eqn{\gamma} represented by basis coefficients `c` on the uniform grid.
#'
#' @param c a \eqn{D} by \eqn{N} matrix (or a vector) of basis coefficients;
#'   each column is an independent warp
#' @param Phi a \eqn{T} by \eqn{D} basis matrix, see [bayes_warp_basis()]
#' @param rep `0` for the shooting vector (v-transform, requires
#'   \eqn{|v| < \pi}), `1` for the clr (h-transform, unconstrained)
#'   representation
#' @param q1 target SRVF, a \eqn{T}-vector shared by all columns or a
#'   \eqn{T} by \eqn{N} matrix
#' @param q2 SRVFs to be warped, a \eqn{T} by \eqn{N} matrix
#' @param var1 noise variance
#' @param grad whether to compute the gradient
#'
#' @return A list with `nll`, `grad`, `SSE`, `gam` and `q2w`.
#'
#' @keywords internal
bayes_warp_geom <- function(c, Phi, rep, q1, q2, var1, grad = TRUE) {
  c <- as.matrix(c)
  q1 <- as.matrix(q1)
  q2 <- as.matrix(q2)
  bayes_warp_geom_c(c, Phi, as.integer(rep), q1, q2, var1, grad)
}

#' One preconditioned (infinite-dimensional) HMC transition
#'
#' Samples \eqn{\pi(c) \propto \exp(-\mathrm{nll}(c)) N(c; 0, \mathrm{diag}(C))}
#' with the splitting integrator of Beskos et al. (2011): the Gaussian prior
#' is integrated exactly by a rotation and the likelihood enters through half
#' kicks \eqn{v \leftarrow v - \epsilon/2\, C \nabla \mathrm{nll}}. The
#' integrator is reversible and volume preserving, so acceptance uses the
#' exact change in \eqn{H = \mathrm{nll}(c) + c'C^{-1}c/2 + v'C^{-1}v/2}.
#'
#' @param c current state
#' @param nll,g negative log-likelihood and gradient of `geom` at `c`
#' @param geom a function of `c` returning a list with `nll` and `grad`
#' @param C prior variances
#' @param epsilon step size (rotation angle)
#' @param L maximum number of leapfrog steps (the number used is uniform on
#'   `1:L`)
#'
#' @return A list with the new state `c`, its `nll` and `g`, whether the
#'   proposal was accepted, the number of leapfrog steps used and the
#'   acceptance probability `alpha` (`0` if the trajectory left the support),
#'   for step-size adaptation.
#'
#' @references Beskos, A., Pinski, F. J., Sanz-Serna, J. M., and Stuart, A. M.
#'   (2011). Hybrid Monte Carlo on Hilbert spaces. Stochastic Processes and
#'   their Applications, 121(10), 2201-2230.
#'
#' @keywords internal
infhmc_step <- function(c, nll, g, geom, C, epsilon, L) {
  C <- as.numeric(C)
  v <- sqrt(C) * stats::rnorm(length(c))
  H0 <- nll + 0.5 * sum(c ^ 2 / C) + 0.5 * sum(v ^ 2 / C)

  nsteps <- ceiling(stats::runif(1) * L)
  q <- c
  gq <- g
  cs <- cos(epsilon)
  sn <- sin(epsilon)
  nll_q <- nll
  for (l in 1:nsteps) {
    v <- v - epsilon / 2 * C * gq
    q_new <- q * cs + v * sn
    v <- v * cs - q * sn
    q <- q_new
    out <- geom(q)
    nll_q <- out$nll
    gq <- out$grad
    if (!is.finite(nll_q)) break
    v <- v - epsilon / 2 * C * gq
  }

  H1 <- nll_q + 0.5 * sum(q ^ 2 / C) + 0.5 * sum(v ^ 2 / C)
  logRatio <- H0 - H1
  alpha <- if (is.finite(logRatio)) min(1, exp(logRatio)) else 0
  accept <- is.finite(logRatio) && (log(stats::runif(1)) < logRatio)
  if (accept) {
    c <- q
    nll <- nll_q
    g <- gq
  }
  list(c = c, nll = nll, g = g, accept = accept, nsteps = nsteps,
       alpha = alpha)
}

#' Preconditioned HMC transitions for independent blocks
#'
#' Runs [infhmc_step()] on \eqn{N} conditionally independent blocks at once
#' (one column of `c` per block), sharing a target of the form
#' \eqn{\prod_i \exp(-\mathrm{nll}_i(c_i)) N(c_i; 0, \mathrm{diag}(C))}. All
#' blocks integrate the same (random) number of leapfrog steps, with their
#' own step sizes, and each is accepted or rejected on its own energy change.
#'
#' @param c a \eqn{D} by \eqn{N} matrix of current states
#' @param nll,g negative log-likelihoods and gradients of `geom` at `c`
#' @param aux extra output of `geom` at `c` (a matrix with \eqn{N} columns,
#'   or `NULL`)
#' @param geom a function of `c` returning a list with `nll`, `grad` and
#'   optionally `aux`, evaluated for all columns at once
#' @param C prior variances
#' @param epsilon step sizes, scalar or of length \eqn{N}
#' @param L maximum number of leapfrog steps
#'
#' @return A list with the new states and their geometry, the per-block
#'   acceptance indicators and the per-block acceptance probabilities `alpha`
#'   (`0` where the trajectory left the support), for step-size adaptation.
#'
#' @keywords internal
infhmc_step_batch <- function(c, nll, g, aux, geom, C, epsilon, L) {
  C <- as.numeric(C)
  D <- nrow(c)
  N <- ncol(c)
  epsilon <- rep_len(as.numeric(epsilon), N)
  v <- matrix(sqrt(C) * stats::rnorm(D * N), D, N)
  H0 <- nll + 0.5 * colSums(c ^ 2 / C) + 0.5 * colSums(v ^ 2 / C)

  nsteps <- ceiling(stats::runif(1) * L)
  q <- c
  gq <- g
  cs <- rep(cos(epsilon), each = D)
  sn <- rep(sin(epsilon), each = D)
  eh <- rep(epsilon / 2, each = D)
  nll_q <- nll
  aux_q <- aux
  failed <- rep(FALSE, N)
  for (l in 1:nsteps) {
    v <- v - eh * C * gq
    q_new <- q * cs + v * sn
    v <- v * cs - q * sn
    q <- q_new
    out <- geom(q)
    nll_q <- out$nll
    gq <- out$grad
    aux_q <- out$aux
    # a block whose trajectory leaves the support (nll = Inf) is rejected, as
    # in infhmc_step(); zero its force so it stays finite
    bad <- !is.finite(nll_q)
    failed <- failed | bad
    if (any(bad)) gq[, bad] <- 0
    v <- v - eh * C * gq
  }

  H1 <- nll_q + 0.5 * colSums(q ^ 2 / C) + 0.5 * colSums(v ^ 2 / C)
  logRatio <- H0 - H1
  ok <- !failed & is.finite(logRatio)
  alpha <- rep(0, N)
  alpha[ok] <- pmin(1, exp(logRatio[ok]))
  accept <- ok & (log(stats::runif(N)) < logRatio)
  if (any(accept)) {
    c[, accept] <- q[, accept]
    nll[accept] <- nll_q[accept]
    g[, accept] <- gq[, accept]
    if (!is.null(aux) && !is.null(aux_q))
      aux[, accept] <- aux_q[, accept]
  }
  list(c = c, nll = nll, g = g, aux = aux, accept = accept, alpha = alpha)
}

#' B-spline basis matrix
#'
#' Evaluates a B-spline basis of `nbasis` functions of order `norder` on
#' `length(x)` equally spaced points spanning the range of `x`. The knot
#' sequence is `nbasis - norder + 2` equally spaced breakpoints with
#' `norder`-fold multiplicity at each end of the range.
#'
#' @param x vector of argument values; only its range and length are used
#' @param nbasis number of basis functions
#' @param norder order of the B-splines, i.e. degree + 1 (default: 4, cubic)
#'
#' @return A `length(x)` by `nbasis` matrix of basis function values.
#'
#' @keywords internal
create_basismatrix <- function(x, nbasis, norder = 4) {
  nbreaks <- nbasis - norder + 2
  if (nbreaks < 2)
    stop("nbasis must be at least norder")
  rangex <- range(x)
  breaks <- seq(rangex[1], rangex[2], length.out = nbreaks)
  knots <- c(rep(rangex[1], norder - 1), breaks, rep(rangex[2], norder - 1))
  xx <- seq(rangex[1], rangex[2], length.out = length(x))
  splines::splineDesign(knots, xx, ord = norder, outer.ok = TRUE)
}

# (T,T) matrix mapping f to the derivative, at the grid, of its interpolating
# cubic spline, as computed in `f_to_srvf()`
spline_deriv_matrix <- function(t) {
  TT <- length(t)
  Dspl <- matrix(0, TT, TT)
  for (k in 1:TT) {
    e <- rep(0, TT)
    e[k] <- 1
    Dspl[, k] <- stats::splinefun(t, e, method = "fmm")(t, deriv = 1)
  }
  Dspl
}

# Inverse of the Karcher mean of warping functions `gam` (T,N) on the uniform
# grid `t`; same estimate as `invertGamma(SqrtMean(gam)$gam_mu)`, vectorized
# over the curves and started from the extrinsic mean
karcher_mean_inverse <- function(gam, t) {
  TT <- nrow(gam)
  N <- ncol(gam)
  dt <- t[2] - t[1]
  w <- rep(dt, TT)
  w[c(1, TT)] <- dt / 2
  gd <- apply(gam, 2, gradient, binsize = dt)
  psi <- sqrt(pmax(gd, 0))
  psi <- sweep(psi, 2, sqrt(colSums(w * psi ^ 2)), "/")
  mu <- rowMeans(psi)
  mu <- mu / sqrt(sum(w * mu ^ 2))
  for (it in 1:100) {
    ip <- pmin(pmax(colSums(w * psi * mu), -1), 1)
    th <- acos(ip)
    fac <- rep(1, N)
    nz <- th > 1e-10
    fac[nz] <- th[nz] / sin(th[nz])
    vbar <- rowMeans(sweep(psi - outer(mu, ip), 2, fac, "*"))
    nv <- sqrt(sum(w * vbar ^ 2))
    if (nv < 1e-8) break
    mu <- cos(nv) * mu + sin(nv) / nv * vbar
    mu <- mu / sqrt(sum(w * mu ^ 2))
  }
  gam_mu <- cumtrapz(t, mu ^ 2)
  gam_mu <- (gam_mu - gam_mu[1]) / (gam_mu[TT] - gam_mu[1])
  as.numeric(invertGamma(gam_mu))
}


#' Dual-averaging step-size adaptation for HMC
#'
#' Nesterov dual averaging as in Hoffman and Gelman (2014), "The No-U-Turn
#' Sampler", Algorithm 5. Call [dual_avg_init()] once to initialize, then
#' this function once per burn-in iteration with the acceptance probability
#' of the transition; at the end of burn-in use `exp(state$log_eps_bar)`, the
#' averaged step size, which is much less noisy than the last iterate. Works
#' elementwise for a vector of step sizes.
#'
#' @param state adaptation state, from [dual_avg_init()]
#' @param alpha acceptance probability (or probabilities) of the last
#'   transition, in \eqn{[0,1]}
#' @param target target acceptance probability
#'
#' @return A list with the step size(s) `epsilon` for the next iteration
#'   (capped at \eqn{\pi/2}) and the updated `state`, whose `log_eps_bar`
#'   holds the averaged log step size.
#'
#' @references Hoffman, M. D. and Gelman, A. (2014). The No-U-Turn sampler:
#'   adaptively setting path lengths in Hamiltonian Monte Carlo. Journal of
#'   Machine Learning Research, 15, 1593-1623.
#'
#' @keywords internal
dual_avg_step <- function(state, alpha, target) {
  t0 <- 10
  gamma <- 0.05
  kappa <- 0.75

  state$n <- state$n + 1
  n <- state$n
  alpha <- as.numeric(alpha)
  state$Hbar <- (1 - 1 / (n + t0)) * state$Hbar + (target - alpha) / (n + t0)
  log_eps <- pmin(state$mu - sqrt(n) / gamma * state$Hbar, log(pi / 2))
  eta <- n ^ (-kappa)
  state$log_eps_bar <- eta * log_eps + (1 - eta) * state$log_eps_bar
  list(epsilon = exp(log_eps), state = state)
}

#' Initialize the dual-averaging step-size adaptation
#'
#' @param eps0 initial step size(s)
#'
#' @return The initial adaptation state for [dual_avg_step()].
#'
#' @keywords internal
dual_avg_init <- function(eps0) {
  list(mu = log(10 * eps0), Hbar = rep(0, length(eps0)),
       log_eps_bar = log(eps0), n = 0)
}
