# Pairwise Bayesian alignment with infinite-dimensional HMC.
#
# Ported from the MATLAB package (pairwise_align_bayes_infHMCbasis.m).

#' Align two functions using hierarchical Bayesian framework with
#' infinite-dimensional HMC
#'
#' This function aligns two functions using a hierarchical Bayesian framework
#' assuming measurement error. It will align `f2` to `f1`. It is based on
#' mapping warping functions to a hypersphere, and a subsequent exponential
#' mapping to a tangent space. In the tangent space, the
#' \eqn{\infty}-dimensional HMC algorithm is used to explore both local and
#' global structure in the posterior distribution.
#'
#' With `warp_rep = "shooting"`, \eqn{v = \Phi c} is a shooting vector at the
#' identity of the Hilbert sphere of SRSFs of warping functions and must
#' satisfy \eqn{|v| < \pi} (draws outside are rejected); with
#' `warp_rep = "clr"`, \eqn{h = \Phi c} is the centered log-derivative,
#' \eqn{\gamma' = \exp(h) / \int \exp(h)}, which is unconstrained. The
#' likelihood gradient with respect to \eqn{c} is exact and \eqn{v} is updated
#' with the preconditioned HMC of Beskos et al. (2011). `timet` is assumed to
#' be uniformly sampled.
#'
#' @param f1 numeric vector of \eqn{M} samples of function 1
#' @param f2 numeric vector of \eqn{M} samples of function 2
#' @param timet time vector of length \eqn{M}
#' @param iter length of the chain
#' @param burnin number of burn-in MCMC iterations; defaults to
#'   `min(5e3, floor(iter / 2))`
#' @param thin thinning of the stored posterior draws
#' @param nchains number of MCMC chains to run
#' @param alpha0,beta0 inverse gamma prior parameters for the measurement
#'   error variances
#' @param alpha,beta inverse gamma prior parameters for the SRVF model
#'   variance
#' @param h the initial HMC step size is `sqrt(h)` (a rotation angle)
#' @param L maximum number of leapfrog steps
#' @param adapt whether to adapt the HMC step size during burn-in
#' @param target_accept target HMC acceptance rate when adapting
#' @param f1propvar,f2propvar proposal variances for `f1` and `f2`
#' @param L1propvar,L2propvar proposal standard deviations of the Gaussian
#'   process length scales of `f1` and `f2`
#' @param vpriorvar prior variance scale of the warping coefficients
#' @param vpriordecay the prior variance of coefficient pair \eqn{j} is
#'   `vpriorvar / j^vpriordecay`
#' @param init.coef initial coefficients of the warping function; length must
#'   be even
#' @param basis basis for the warping coefficients, either `"fourier"` or
#'   `"legendre"`
#' @param warp_rep warping function representation, either `"shooting"` (the
#'   v-transform) or `"clr"` (the h-transform)
#' @param init where the first chain starts: `"dp"` at the elastic (dynamic
#'   programming) alignment of the smoothed functions, or `"coef"` at
#'   `init.coef`. Further chains always start at random coefficients. A chain
#'   started at the identity can settle in a poor local mode
#' @param sampfreq subsampling frequency used when fitting the initial
#'   Gaussian process models to `f1` and `f2`
#' @param extrainfo whether additional information is returned
#' @param verbose whether to show a progress bar
#'
#' @return Returns a list containing
#' \item{f1}{posterior draws of \eqn{f_1}}
#' \item{f2}{posterior draws of \eqn{f_2}}
#' \item{gamma}{posterior mean warping function}
#' \item{v_coef}{posterior draws of the warping coefficients}
#' \item{psi}{posterior mean psi function}
#' \item{sigma}{posterior draws of the SRVF model standard deviation}
#' \item{sigma1,sigma2}{posterior draws of the measurement error variances}
#' \item{s1,s2}{posterior draws of the Gaussian process variances}
#' \item{L1,L2}{posterior draws of the Gaussian process length scales}
#' \item{f2_warped_mu}{`f2` warped by the posterior mean warping function}
#' \item{epsilon}{HMC step size at each iteration}
#' \item{theta_accept}{acceptance of the warping samples (if `extrainfo`)}
#' \item{f1_accept,f2_accept}{acceptance of the `f1`, `f2` samples (if
#'   `extrainfo`)}
#' \item{L1_accept,L2_accept}{acceptance of the length scales (if
#'   `extrainfo`)}
#' \item{gamma_mat}{matrix of posterior draws of gamma (if `extrainfo`)}
#' \item{gamma_stats}{pointwise 95% credible band of gamma (if `extrainfo`)}
#' \item{SSE}{SSE of the level 2 likelihood (if `extrainfo`)}
#' \item{logl}{log-likelihood trace (if `extrainfo`)}
#' \item{f2_warped}{posterior draws of the aligned `f2` (if `extrainfo`)}
#' \item{phasedist}{posterior draws of the phase distance (if `extrainfo`)}
#' \item{ampdist}{posterior draws of the amplitude distance (if `extrainfo`)}
#' \item{posterior_gamma_modes}{modes of the posterior of gamma, found by
#'   clustering the chains (if `nchains > 1`)}
#' \item{posterior_gamma_modes_cr}{credible regions of the modes (if
#'   `nchains > 1`)}
#'
#' @keywords bayesian
#' @concept srsf alignment
#'
#' @references Tucker, J. D., Shand, L., and Chowdhary, K. (2021). Multimodal
#'   Bayesian registration of noisy functions using Hamiltonian Monte Carlo.
#'   Computational Statistics and Data Analysis, 163, 107298.
#' @references Beskos, A., Pinski, F. J., Sanz-Serna, J. M., and Stuart, A. M.
#'   (2011). Hybrid Monte Carlo on Hilbert spaces. Stochastic Processes and
#'   their Applications, 121(10), 2201-2230.
#'
#' @export
#' @examples
#' \dontrun{
#'   # This is an MCMC algorithm and takes a long time to run
#'   out <- pair_align_functions_infHMC(
#'     f1 = simu_data$f[, 1],
#'     f2 = simu_data$f[, 21],
#'     timet = simu_data$time,
#'     iter = 1000
#'   )
#'   mean(out$theta_accept)
#' }
pair_align_functions_infHMC <- function(f1,
                                        f2,
                                        timet,
                                        iter = 1e4,
                                        burnin = NULL,
                                        thin = 50,
                                        nchains = 1,
                                        alpha0 = 0.1,
                                        beta0 = 0.1,
                                        alpha = 1,
                                        beta = 1,
                                        h = 0.01,
                                        L = 4L,
                                        adapt = TRUE,
                                        target_accept = 0.65,
                                        f1propvar = 1e-4,
                                        f2propvar = 1e-4,
                                        L1propvar = 0.3,
                                        L2propvar = 0.3,
                                        vpriorvar = 1,
                                        vpriordecay = 1,
                                        init.coef = rep(0, 20),
                                        basis = c("fourier", "legendre"),
                                        warp_rep = c("shooting", "clr"),
                                        init = c("dp", "coef"),
                                        sampfreq = 5L,
                                        extrainfo = FALSE,
                                        verbose = TRUE) {
  basis <- rlang::arg_match(basis)
  warp_rep <- rlang::arg_match(warp_rep)
  init <- rlang::arg_match(init)

  if (length(f1) != length(f2))
    stop("Length of f1 and f2 must be equal")
  if (length(f1) != length(timet))
    stop("Length of f1 and timet must be equal")
  if (length(init.coef) %% 2 != 0)
    stop("Length of init.coef must be even")
  if (is.null(burnin))
    burnin <- min(5e3, floor(iter / 2))
  if (burnin < 1 || burnin >= iter)
    stop("burnin must be between 1 and iter - 1")

  ncoef <- length(init.coef)

  # random start points for more than one chain
  random_starts <- matrix(0, ncoef, nchains)
  if (nchains > 1) {
    for (i in 1:nchains) {
      randcoef <- stats::runif(ncoef, -1, 1)
      if (warp_rep == "shooting" && sqrt(sum(randcoef ^ 2)) >= 0.9 * pi)
        randcoef <- randcoef * (0.9 * pi / sqrt(sum(randcoef ^ 2))) *
          stats::runif(1)
      random_starts[, i] <- randcoef
    }
  }
  # chain 1 starts at init.coef, or at the DP warp when init = "dp": chains
  # started at the identity can settle in a poor local mode
  random_starts[, 1] <- init.coef

  opts <- list(
    iter = iter, burnin = burnin, alpha0 = alpha0, beta0 = beta0,
    alpha = alpha, beta = beta, h = h, L = L, adapt = adapt,
    target_accept = target_accept, f1propvar = f1propvar,
    f2propvar = f2propvar, L1propvar = L1propvar, L2propvar = L2propvar,
    vpriorvar = vpriorvar, vpriordecay = vpriordecay, basis = basis,
    warp_rep = warp_rep, init = init, sampfreq = sampfreq,
    extrainfo = extrainfo,
    verbose = verbose && nchains == 1
  )

  chains <- vector("list", nchains)
  for (i in 1:nchains) {
    opts$initcoef <- random_starts[, i]
    opts$dp_start <- (i == 1) && init == "dp"
    chains[[i]] <- run_mcmc_infHMC(f1, f2, timet, opts)
  }

  out <- combine_infHMC_chains(chains, nchains, extrainfo)

  # finding modes
  if (nchains > 1) {
    M <- nrow(out$gamma)
    Dx <- matrix(0, nchains, nchains)
    time1 <- seq(0, 1, length.out = M)
    binsize <- mean(diff(time1))
    for (i in 1:(nchains - 1)) {
      for (j in (i + 1):nchains) {
        psi1 <- sqrt(pmax(gradient(out$gamma[, i], binsize), 0))
        psi2 <- sqrt(pmax(gradient(out$gamma[, j], binsize), 0))
        Dx[i, j] <- acos(min(max(trapz(time1, psi1 * psi2), -1), 1))
      }
    }
    Dx <- Dx + t(Dx)

    Z <- stats::hclust(stats::as.dist(Dx), method = "complete")
    cutoff <- stats::median(Dx)
    cl <- stats::cutree(Z, h = cutoff)
    nmodes <- length(unique(cl))

    out$posterior_gamma_modes <- matrix(0, M, nmodes)
    out$posterior_gamma_modes_cr <- array(0, dim = c(M, 2, nmodes))
    for (i in 1:nmodes) {
      idx <- which(cl == sort(unique(cl))[i])
      tmp <- do.call(cbind, lapply(idx, function(k) chains[[k]]$gamma_mat))
      out$posterior_gamma_modes[, i] <- SqrtMean(tmp)$gam_mu
      out$posterior_gamma_modes_cr[, , i] <- t(statsFun_gam(tmp))
    }
  }

  # thinning
  keep <- seq(1, nrow(out$f1), by = thin)
  out$f1 <- out$f1[keep, , drop = FALSE]
  out$f2 <- out$f2[keep, , drop = FALSE]
  out$v_coef <- out$v_coef[keep, , drop = FALSE]
  for (nm in c("sigma", "sigma1", "sigma2", "s1", "s2", "L1", "L2"))
    out[[nm]] <- out[[nm]][keep]

  if (extrainfo) {
    keepa <- seq(1, length(out$theta_accept), by = thin)
    for (nm in c("theta_accept", "f1_accept", "f2_accept", "L1_accept",
                 "L2_accept"))
      out[[nm]] <- out[[nm]][keepa]
    out$gamma_mat <- out$gamma_mat[, keep, drop = FALSE]
    keeps <- seq(1, length(out$SSE), by = thin)
    out$SSE <- out$SSE[keeps]
    out$logl <- out$logl[keeps]
    out$f2_warped <- out$f2_warped[keep, , drop = FALSE]
    out$phasedist <- out$phasedist[keep]
    out$ampdist <- out$ampdist[keep]
  }

  out
}

statsFun_gam <- function(mat)
  apply(mat, 1, stats::quantile, probs = c(0.025, 0.975))

combine_infHMC_chains <- function(chains, nchains, extrainfo) {
  Nsamples <- nrow(chains[[1]]$f1)
  M <- ncol(chains[[1]]$f1)
  out <- list()
  out$f1 <- do.call(rbind, lapply(chains, `[[`, "f1"))
  out$f2 <- do.call(rbind, lapply(chains, `[[`, "f2"))
  out$gamma <- do.call(cbind, lapply(chains, `[[`, "gamma"))
  out$v_coef <- do.call(rbind, lapply(chains, `[[`, "v_coef"))
  out$psi <- do.call(cbind, lapply(chains, `[[`, "psi"))
  for (nm in c("sigma", "sigma1", "sigma2", "s1", "s2", "L1", "L2"))
    out[[nm]] <- unlist(lapply(chains, `[[`, nm))
  out$f2_warped_mu <- do.call(cbind, lapply(chains, `[[`, "f2_warped_mu"))
  out$epsilon <- do.call(rbind, lapply(chains, `[[`, "epsilon"))

  if (extrainfo) {
    for (nm in c("theta_accept", "f1_accept", "f2_accept", "L1_accept",
                 "L2_accept", "SSE", "logl", "phasedist", "ampdist"))
      out[[nm]] <- unlist(lapply(chains, `[[`, nm))
    out$gamma_mat <- do.call(cbind, lapply(chains, `[[`, "gamma_mat"))
    out$gamma_stats <- t(statsFun_gam(out$gamma_mat))
    out$f2_warped <- do.call(rbind, lapply(chains, `[[`, "f2_warped"))
  }

  if (nchains == 1) {
    out$gamma <- as.numeric(out$gamma)
    out$psi <- as.numeric(out$psi)
    out$f2_warped_mu <- as.numeric(out$f2_warped_mu)
    out$epsilon <- as.numeric(out$epsilon)
  }

  out
}

# squared-exponential correlation function
exp2corr <- function(sigma2, phi, ds) sigma2 * exp(-ds ^ 2 / (2 * phi ^ 2))

# GP prior N(0, s * (corr(L) + 0.1 I)) stored through its Cholesky factor
# (Dsq: squared pairwise distances of the time points)
gp_prior <- function(s, L, Dsq, Rcorr = NULL) {
  n <- nrow(Dsq)
  if (is.null(Rcorr)) {
    K <- exp(-Dsq / (2 * L ^ 2))
    diag(K) <- diag(K) + 0.1
    Rcorr <- chol(K)
  }
  list(Rcorr = Rcorr, s = s,
       logdet = 2 * sum(log(diag(Rcorr))) + n * log(s))
}

quad_form <- function(R, f) {
  z <- backsolve(R, as.numeric(f), transpose = TRUE)
  sum(z ^ 2)
}

chol_jitter <- function(K) {
  n <- nrow(K)
  jit <- 0
  R <- try(chol(K), silent = TRUE)
  while (inherits(R, "try-error")) {
    jit <- max(2 * jit, 1e-10 * mean(diag(K)))
    R <- try(chol(K + jit * diag(n)), silent = TRUE)
  }
  R
}

# Metropolis update of a latent function f (f1 or f2) given the data y, its
# GP prior and the SRVF alignment likelihood. f | y, prior is Gaussian with
# precision K^{-1} + I/sigmaf, so the target is that density times the
# alignment likelihood exp(-SSE/(2 sigma)).
update_f_infHMC <- function(f_curr, q_curr, SSE_curr, y, sse_fun, gp, Rprop,
                            to_srvf, sigma, sigmaf) {
  f_prop <- f_curr + as.numeric(stats::rnorm(length(f_curr)) %*% Rprop)
  q_prop <- to_srvf(f_prop)
  SSE_prop <- sse_fun(q_prop)
  lp_curr <- -SSE_curr / (2 * sigma) -
    (quad_form(gp$Rcorr, f_curr) / gp$s + sum((y - f_curr) ^ 2) / sigmaf) / 2
  lp_prop <- -SSE_prop / (2 * sigma) -
    (quad_form(gp$Rcorr, f_prop) / gp$s + sum((y - f_prop) ^ 2) / sigmaf) / 2
  accept <- log(stats::runif(1)) < (lp_prop - lp_curr)
  if (is.na(accept)) accept <- FALSE
  if (accept) {
    f_curr <- f_prop
    q_curr <- q_prop
    SSE_curr <- SSE_prop
  }
  list(f = f_curr, q = q_curr, SSE = SSE_curr, accept = accept)
}

# Metropolis-Hastings update of a GP length scale with a truncated normal
# random-walk proposal (flat prior on L > 0)
update_length_infHMC <- function(f, gp, s, L_curr, propsd, Dsq) {
  # truncated normal N(L_curr, propsd^2) on (0, inf) by rejection (exact)
  L_prop <- L_curr + propsd * stats::rnorm(1)
  while (L_prop <= 0)
    L_prop <- L_curr + propsd * stats::rnorm(1)
  gp_prop <- gp_prior(s, L_prop, Dsq)
  lp_curr <- -gp$logdet / 2 - quad_form(gp$Rcorr, f) / (2 * s)
  lp_prop <- -gp_prop$logdet / 2 - quad_form(gp_prop$Rcorr, f) / (2 * s)
  # Hastings correction for the truncation: q(L'|L) = phi(.)/Phi(L/propsd)
  log_hastings <- log(2 * stats::pnorm(L_curr / propsd)) -
    log(2 * stats::pnorm(L_prop / propsd))
  accept <- log(stats::runif(1)) < (lp_prop - lp_curr + log_hastings)
  if (is.na(accept)) accept <- FALSE
  if (accept) {
    L_curr <- L_prop
    gp <- gp_prop
  }
  list(L = L_curr, gp = gp, accept = accept)
}

# Maximum (marginal) likelihood fit of a squared-exponential Gaussian process
# regression with a constant mean and Gaussian noise, the R counterpart of
# MATLAB's fitrgp with its default options. Returns the kernel length scale,
# the kernel variance and the posterior mean at `xout`.
fit_gp_sqexp <- function(x, y, xout = x) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(y)
  Dsq <- as.matrix(stats::dist(x)) ^ 2
  sdy <- stats::sd(y)
  if (sdy <= 0) sdy <- 1
  p0 <- log(c(max(stats::sd(x), .Machine$double.eps), sdy / sqrt(2),
              sdy / sqrt(2)))

  nll <- function(p) {
    L <- exp(p[1]); sf <- exp(p[2]); sn <- exp(p[3])
    K <- sf ^ 2 * exp(-Dsq / (2 * L ^ 2))
    diag(K) <- diag(K) + sn ^ 2 + 1e-10
    R <- try(chol(K), silent = TRUE)
    if (inherits(R, "try-error")) return(1e10)
    one <- rep(1, n)
    Ky <- backsolve(R, backsolve(R, y, transpose = TRUE))
    K1 <- backsolve(R, backsolve(R, one, transpose = TRUE))
    mu <- sum(one * Ky) / sum(one * K1)
    z <- backsolve(R, y - mu, transpose = TRUE)
    val <- sum(log(diag(R))) + sum(z ^ 2) / 2
    if (!is.finite(val)) 1e10 else val
  }

  fit <- try(stats::optim(p0, nll, method = "BFGS",
                          control = list(maxit = 200)), silent = TRUE)
  p <- if (inherits(fit, "try-error")) p0 else fit$par
  L <- exp(p[1]); sf <- exp(p[2]); sn <- exp(p[3])

  K <- sf ^ 2 * exp(-Dsq / (2 * L ^ 2))
  diag(K) <- diag(K) + sn ^ 2 + 1e-10
  R <- chol(K)
  one <- rep(1, n)
  Ky <- backsolve(R, backsolve(R, y, transpose = TRUE))
  K1 <- backsolve(R, backsolve(R, one, transpose = TRUE))
  mu <- sum(one * Ky) / sum(one * K1)
  wts <- backsolve(R, backsolve(R, y - mu, transpose = TRUE))
  Ks <- sf ^ 2 * exp(-outer(as.numeric(xout), x, "-") ^ 2 / (2 * L ^ 2))
  list(L = L, s = sf ^ 2, sigma_n = sn, pred = as.numeric(mu + Ks %*% wts))
}

# main mcmc function
run_mcmc_infHMC <- function(y1i, y2i, timet, opts) {
  SIG_GAM <- 13
  iter <- opts$iter
  y1 <- as.numeric(y1i)
  y2 <- as.numeric(y2i)
  timet <- as.numeric(timet)
  TT <- length(timet)

  # normalize time to [0,1] (uniform sampling assumed)
  t01 <- (timet - min(timet)) / (max(timet) - min(timet))

  # warping function parameter settings
  valid_index <- opts$burnin:iter
  ncoef <- length(opts$initcoef)
  rep <- as.numeric(opts$warp_rep == "clr")
  Phi <- bayes_warp_basis(t01, ncoef, opts$basis)
  C <- opts$vpriorvar / rep(1:(ncoef / 2), each = 2) ^ opts$vpriordecay
  epsilon <- sqrt(opts$h)
  da <- dual_avg_init(epsilon)
  L <- opts$L

  # f1, f2 GP priors, initialized from a GP fit
  Dmat <- as.matrix(stats::dist(timet))
  Dsq <- Dmat ^ 2

  # SRVF operator: the derivative of the interpolating cubic spline used by
  # f_to_srvf() is linear in f, so build its (T,T) matrix once; to_srvf(f) is
  # identical to f_to_srvf(f, t01) without fitting a spline every update
  Dspl <- spline_deriv_matrix(t01)
  eps <- .Machine$double.eps
  to_srvf <- function(f) {
    fy <- as.numeric(Dspl %*% as.numeric(f))
    fy / sqrt(abs(fy) + eps)
  }

  sub <- seq(1, TT, by = opts$sampfreq)
  model1 <- fit_gp_sqexp(timet[sub], y1[sub], timet)
  model2 <- fit_gp_sqexp(timet[sub], y2[sub], timet)

  s1_curr <- model1$s
  s2_curr <- model2$s
  L1_curr <- model1$L
  L2_curr <- model2$L
  gp1 <- gp_prior(s1_curr, L1_curr, Dsq)
  gp2 <- gp_prior(s2_curr, L2_curr, Dsq)

  # random-walk proposals for f1, f2 (fixed covariance, factor once)
  R_f1prop <- chol_jitter(exp2corr(opts$f1propvar, L1_curr, Dmat))
  R_f2prop <- chol_jitter(exp2corr(opts$f2propvar, L2_curr, Dmat))

  # result objects
  res_v_coef <- matrix(0, iter, ncoef)
  res_sigma <- res_sigma1 <- res_sigma2 <- numeric(iter)
  res_f1 <- res_f2 <- matrix(0, iter, TT)
  res_s1 <- res_s2 <- res_L1 <- res_L2 <- numeric(iter)
  res_logl <- res_SSE <- res_epsilon <- numeric(iter)
  res_theta_accept <- res_f1_accept <- res_f2_accept <- numeric(iter)
  res_L1_accept <- res_L2_accept <- numeric(iter)

  # init
  v_coef_curr <- as.numeric(opts$initcoef)
  sigma_curr <- 1
  sigma1_curr <- 0.01
  sigma2_curr <- 0.01
  f1_curr <- model1$pred
  f2_curr <- model2$pred
  q1_curr <- to_srvf(f1_curr)
  q2_curr <- to_srvf(f2_curr)
  if (isTRUE(opts$dp_start)) {
    # start at the elastic (DP) alignment of the GP-smoothed functions
    gam_dp <- optimum.reparam(q1_curr, t01, q2_curr, t01)
    v_coef_curr <- as.numeric(bayes_gam_to_coef(gam_dp, Phi, rep))
  }

  geo <- bayes_warp_geom(v_coef_curr, Phi, rep, q1_curr, q2_curr, sigma_curr)
  nll <- geo$nll[1]
  SSE_curr <- geo$SSE[1]
  if (!is.finite(nll))
    stop("Initial warping coefficients are outside the injectivity radius (|v| >= pi)")

  res_v_coef[1, ] <- v_coef_curr
  res_f1[1, ] <- f1_curr
  res_f2[1, ] <- f2_curr
  res_sigma[1] <- sigma_curr
  res_sigma1[1] <- sigma1_curr
  res_sigma2[1] <- sigma2_curr
  res_s1[1] <- s1_curr
  res_s2[1] <- s2_curr
  res_L1[1] <- L1_curr
  res_L2[1] <- L2_curr
  res_SSE[1] <- SSE_curr
  res_logl[1] <- -nll
  res_epsilon[1] <- epsilon

  pb <- NULL
  if (opts$verbose)
    pb <- cli::cli_progress_bar("Running MCMC", total = iter - 1,
                                .envir = environment())

  for (m in 2:iter) {
    # update f1 | v, f2
    upd <- update_f_infHMC(
      f1_curr, q1_curr, SSE_curr, y1,
      function(q) bayes_warp_geom(v_coef_curr, Phi, rep, q, q2_curr,
                                  sigma_curr, grad = FALSE)$SSE[1],
      gp1, R_f1prop, to_srvf, sigma_curr, sigma1_curr)
    f1_curr <- upd$f; q1_curr <- upd$q; SSE_curr <- upd$SSE
    f1_accept <- upd$accept

    # update f2 | v, f1
    upd <- update_f_infHMC(
      f2_curr, q2_curr, SSE_curr, y2,
      function(q) bayes_warp_geom(v_coef_curr, Phi, rep, q1_curr, q,
                                  sigma_curr, grad = FALSE)$SSE[1],
      gp2, R_f2prop, to_srvf, sigma_curr, sigma2_curr)
    f2_curr <- upd$f; q2_curr <- upd$q
    f2_accept <- upd$accept

    # update v by infinite-dimensional HMC; geometry evaluated at the current
    # f1, f2 and sigma
    geom <- function(cc) {
      o <- bayes_warp_geom(cc, Phi, rep, q1_curr, q2_curr, sigma_curr)
      list(nll = o$nll[1], grad = as.numeric(o$grad))
    }
    g0 <- geom(v_coef_curr)
    hmc <- infhmc_step(v_coef_curr, g0$nll, g0$grad, geom, C, epsilon, L)
    v_coef_curr <- hmc$c
    theta_accept <- hmc$accept
    SSE_curr <- 2 * sigma_curr * (hmc$nll - TT / 2 * log(sigma_curr))

    # adapt step size during burn-in (dual averaging); keep the averaged
    # step size afterwards
    if (opts$adapt && m <= opts$burnin) {
      da_out <- dual_avg_step(da, hmc$alpha, opts$target_accept)
      epsilon <- da_out$epsilon
      da <- da_out$state
      if (m == opts$burnin)
        epsilon <- exp(da$log_eps_bar)
    }

    # update sigma^2 (SRVF model variance)
    newshape <- TT / 2 + opts$alpha
    newscale <- SSE_curr / 2 + opts$beta
    sigma_curr <- 1 / stats::rgamma(1, shape = newshape, scale = 1 / newscale)

    # update sigma1^2, sigma2^2 (measurement error variances)
    newshape <- TT / 2 + opts$alpha0
    sigma1_curr <- 1 / stats::rgamma(
      1, shape = newshape,
      scale = 1 / (sum((y1 - f1_curr) ^ 2) / 2 + opts$beta0))
    sigma2_curr <- 1 / stats::rgamma(
      1, shape = newshape,
      scale = 1 / (sum((y2 - f2_curr) ^ 2) / 2 + opts$beta0))

    # update GP variances s1^2, s2^2
    s1_curr <- 1 / stats::rgamma(
      1, shape = newshape,
      scale = 1 / (quad_form(gp1$Rcorr, f1_curr) / 2 + opts$beta0))
    s2_curr <- 1 / stats::rgamma(
      1, shape = newshape,
      scale = 1 / (quad_form(gp2$Rcorr, f2_curr) / 2 + opts$beta0))
    gp1 <- gp_prior(s1_curr, L1_curr, Dsq, gp1$Rcorr)
    gp2 <- gp_prior(s2_curr, L2_curr, Dsq, gp2$Rcorr)

    # update GP length scales L1, L2
    ul <- update_length_infHMC(f1_curr, gp1, s1_curr, L1_curr, opts$L1propvar,
                               Dsq)
    L1_curr <- ul$L; gp1 <- ul$gp; L1_accept <- ul$accept
    ul <- update_length_infHMC(f2_curr, gp2, s2_curr, L2_curr, opts$L2propvar,
                               Dsq)
    L2_curr <- ul$L; gp2 <- ul$gp; L2_accept <- ul$accept

    # save update to results
    res_v_coef[m, ] <- v_coef_curr
    res_f1[m, ] <- f1_curr
    res_f2[m, ] <- f2_curr
    res_sigma[m] <- sigma_curr
    res_sigma1[m] <- sigma1_curr
    res_sigma2[m] <- sigma2_curr
    res_s1[m] <- s1_curr
    res_s2[m] <- s2_curr
    res_L1[m] <- L1_curr
    res_L2[m] <- L2_curr
    res_SSE[m] <- SSE_curr
    res_logl[m] <- -(SSE_curr / (2 * sigma_curr) + TT / 2 * log(sigma_curr))
    res_epsilon[m] <- epsilon
    res_theta_accept[m] <- theta_accept
    res_f1_accept[m] <- f1_accept
    res_f2_accept[m] <- f2_accept
    res_L1_accept[m] <- L1_accept
    res_L2_accept[m] <- L2_accept

    if (!is.null(pb)) cli::cli_progress_update(id = pb)
  }
  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  # posterior warping functions
  cg <- bayes_coef_to_gam(t(res_v_coef[valid_index, , drop = FALSE]), Phi, rep)
  gamma_mat <- round(cg$gam, SIG_GAM)
  psi_mat <- cg$psi

  # posterior mean of psi (extrinsic mean, projected back to the sphere)
  psi_mean <- rowMeans(psi_mat)
  psi_mean <- psi_mean / sqrt(trapz(t01, psi_mean ^ 2))
  gam_mean <- c(0, cumsum(diff(t01) * (psi_mean[-TT] ^ 2 + psi_mean[-1] ^ 2) / 2))
  gam_mean <- norm_gam(gam_mean)

  out <- list()
  out$f1 <- res_f1[valid_index, , drop = FALSE]
  out$f2 <- res_f2[valid_index, , drop = FALSE]
  out$gamma <- gam_mean
  out$v_coef <- res_v_coef[valid_index, , drop = FALSE]
  out$psi <- psi_mean
  out$sigma <- res_sigma[valid_index]
  out$sigma1 <- res_sigma1[valid_index]
  out$sigma2 <- res_sigma2[valid_index]
  out$s1 <- res_s1[valid_index]
  out$s2 <- res_s2[valid_index]
  out$L1 <- res_L1[valid_index]
  out$L2 <- res_L2[valid_index]
  out$epsilon <- res_epsilon
  out$f2_warped_mu <- warp_f_gamma(colMeans(out$f2), t01, gam_mean)

  if (opts$extrainfo) {
    nv <- length(valid_index)
    Dx <- Re(acos(pmin(pmax(
      apply(psi_mat, 2, function(p) trapz(t01, p)), -1), 1)))
    Dy <- numeric(nv)
    f2_warped <- matrix(0, nv, TT)
    for (ii in 1:nv) {
      q1s <- to_srvf(out$f1[ii, ])
      q2s <- to_srvf(out$f2[ii, ])
      q2warp <- warp_q_gamma(q2s, t01, gamma_mat[, ii])
      Dy[ii] <- sqrt(trapz(t01, (q1s - q2warp) ^ 2))
      f2_warped[ii, ] <- warp_f_gamma(out$f2[ii, ], t01, gamma_mat[, ii])
    }
    out$theta_accept <- res_theta_accept[-1]
    out$f1_accept <- res_f1_accept[-1]
    out$f2_accept <- res_f2_accept[-1]
    out$L1_accept <- res_L1_accept[-1]
    out$L2_accept <- res_L2_accept[-1]
    out$gamma_mat <- gamma_mat
    out$gamma_stats <- t(statsFun_gam(gamma_mat))
    out$SSE <- res_SSE
    out$logl <- res_logl
    out$phasedist <- Dx
    out$ampdist <- Dy
    out$f2_warped <- f2_warped
  }

  out
}
