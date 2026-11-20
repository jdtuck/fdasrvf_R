# Group-wise Bayesian alignment with infinite-dimensional HMC.
#
# Ported from the MATLAB package (fdawarp.time_warping_bayes_infHMC).

#' Bayesian group-wise alignment using infinite-dimensional HMC
#'
#' Hierarchical Bayesian model in SRVF space. For \eqn{i = 1, \dots, N}
#' \deqn{q^*(t_k) = (q_i \circ \gamma_i)(t_k) \sqrt{\gamma_i'(t_k)} +
#' e_{ik}, \quad e_{ik} \sim N(0, \sigma^2)} with template
#' \eqn{q^* = B b} (B-spline basis, \eqn{b \sim N(0, \mathrm{qpriorvar}\, I)}),
#' warping coefficients \eqn{c_i \sim N(0, \mathrm{diag}(C))},
#' \eqn{C_j = \mathrm{vpriorvar} / j^{\mathrm{vpriordecay}}}, and
#' \eqn{\sigma^2 \sim \mathrm{IG}(\alpha_0, \beta_0)}. The warping functions
#' use either the shooting vector (v-transform, \eqn{|v| < \pi}) or the
#' centered log-derivative (h-transform, unconstrained) representation.
#'
#' Each sweep updates every \eqn{c_i} by preconditioned
#' (infinite-dimensional) HMC with the exact likelihood gradient, then draws
#' \eqn{b} and \eqn{\sigma^2} from their conjugate full conditionals. Step
#' sizes are adapted per curve during burn-in. The amplitude-phase split is
#' identified only up to a common warp (a ridge that per-curve updates move
#' along slowly), so the centering constraint is imposed by recentering the
#' warps every `ncenter` iterations so that their Karcher mean is the identity
#' (the template is warped accordingly). The chain is initialized from an
#' elastic (DP) Karcher mean alignment unless `init = "identity"`.
#'
#' @param f matrix (\eqn{M} x \eqn{N}) of \eqn{N} functions with \eqn{M}
#'   samples
#' @param time vector of length \eqn{M} describing the sample points
#' @param iter number of iterations
#' @param burnin number of burn-in iterations; defaults to `floor(iter / 2)`
#' @param thin thinning of the stored posterior draws
#' @param alpha0,beta0 inverse gamma prior parameters for \eqn{\sigma^2}
#' @param h the initial HMC step size is `sqrt(h)`
#' @param L maximum number of leapfrog steps
#' @param adapt whether to adapt the step sizes during burn-in
#' @param target_accept target acceptance rate when adapting
#' @param ncoef number of warping basis coefficients
#' @param basis warping basis, either `"fourier"` or `"legendre"`
#' @param warp_rep warping function representation, either `"shooting"` (the
#'   v-transform) or `"clr"` (the h-transform)
#' @param vpriorvar prior variance scale for the warping coefficients
#' @param vpriordecay the prior variance decays as `1 / j^vpriordecay`
#' @param nbasis number of B-spline basis functions for \eqn{q^*}
#' @param qpriorvar prior variance of the template coefficients
#' @param init either `"dp"` (initialize from an elastic Karcher mean
#'   alignment) or `"identity"`
#' @param ncenter recenter the warping functions every `ncenter` iterations
#'   (`0` to never recenter)
#' @param center_all whether to also recenter after burn-in (`FALSE`:
#'   during burn-in only)
#' @param extrainfo whether to return posterior draws and summaries
#' @param verbose whether to show a progress bar
#'
#' @return An object of class `fdawarp` which is a list with the same
#'   components as [time_warping()], where `warping_functions` holds the
#'   posterior (extrinsic) mean warping function of each curve, together with
#'   a `mcmc` component containing
#' \item{v_coef}{array (\eqn{D} x \eqn{N} x `nsamp`) of posterior draws of
#'   the warping coefficients}
#' \item{sigma}{posterior draws of \eqn{\sigma} (all iterations)}
#' \item{accept}{matrix (\eqn{N} x `iter - 1`) of HMC acceptance indicators}
#' \item{epsilon}{matrix (\eqn{N} x `iter`) of step sizes}
#' \item{logl}{log-likelihood trace}
#' \item{burnin, iter}{number of burn-in and total iterations}
#' \item{gamma_mat}{array (\eqn{M} x `nsamp` x \eqn{N}) of posterior gammas
#'   (if `extrainfo`)}
#' \item{gamma_stats}{array (\eqn{M} x 2 x \eqn{N}) of pointwise 95% bands of
#'   the gammas (if `extrainfo`)}
#' \item{q_star_mat}{matrix (\eqn{M} x `nsamp`) of posterior templates (if
#'   `extrainfo`)}
#' \item{q_star_stats}{pointwise 95% bands of the template (if `extrainfo`)}
#' \item{fmean_mat}{matrix (\eqn{M} x `nsamp`) of posterior mean functions
#'   (if `extrainfo`)}
#' \item{fmean_stats}{pointwise 95% bands of the mean function (if
#'   `extrainfo`)}
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
#'   out <- time_warping_bayes_infHMC(simu_data$f, simu_data$time, iter = 600)
#' }
time_warping_bayes_infHMC <- function(f,
                                      time,
                                      iter = 2000L,
                                      burnin = NULL,
                                      thin = 1L,
                                      alpha0 = 0.1,
                                      beta0 = 0.1,
                                      h = 0.01,
                                      L = 10L,
                                      adapt = TRUE,
                                      target_accept = 0.65,
                                      ncoef = 10L,
                                      basis = c("fourier", "legendre"),
                                      warp_rep = c("shooting", "clr"),
                                      vpriorvar = 1,
                                      vpriordecay = 1,
                                      nbasis = 40L,
                                      qpriorvar = 1e3,
                                      init = c("dp", "identity"),
                                      ncenter = 20L,
                                      center_all = TRUE,
                                      extrainfo = TRUE,
                                      verbose = TRUE) {
  basis <- rlang::arg_match(basis)
  warp_rep <- rlang::arg_match(warp_rep)
  init <- rlang::arg_match(init)

  f <- as.matrix(f)
  time <- as.numeric(time)
  M <- nrow(f)
  N <- ncol(f)
  if (length(time) != M)
    stop("Length of time must equal nrow(f)")
  if (is.null(burnin))
    burnin <- floor(iter / 2)
  if (burnin < 1 || burnin >= iter)
    stop("burnin must be between 1 and iter - 1")
  if (basis == "fourier" && ncoef %% 2 != 0)
    stop("ncoef must be even for the fourier basis")
  rep <- as.numeric(warp_rep == "clr")

  # work on a uniform grid on [0,1]
  tn <- (time - time[1]) / (time[M] - time[1])
  t01 <- seq(0, 1, length.out = M)
  f0 <- apply(f, 2, function(y) stats::approx(tn, y, xout = t01)$y)
  q <- f_to_srvf(f0, t01)

  Phi <- bayes_warp_basis(t01, ncoef, basis)
  C <- vpriorvar / ceiling((1:ncoef) / 2) ^ vpriordecay
  B <- create_basismatrix(t01, nbasis, 4)
  BtB <- crossprod(B)
  Bpinv <- qr.solve(B, diag(M))

  # initialization
  if (init == "dp") {
    tmp <- if (verbose)
      time_warping(f0, t01, lambda = 0, max_iter = 1L)
    else
      suppressMessages(time_warping(f0, t01, lambda = 0, max_iter = 10L))
    c_curr <- bayes_gam_to_coef(tmp$warping_functions, Phi, rep)
    gam0 <- bayes_coef_to_gam(c_curr, Phi, rep)$gam
    qtil <- matrix(0, M, N)
    for (ii in 1:N)
      qtil[, ii] <- warp_q_gamma(q[, ii], t01, gam0[, ii])
    qstar <- rowMeans(qtil)
  } else {
    c_curr <- matrix(0, ncoef, N)
    dqq <- sqrt(colSums((q - rowMeans(q)) ^ 2))
    qstar <- q[, which.min(dqq)]
  }
  b_curr <- as.numeric(Bpinv %*% qstar)
  qstar <- as.numeric(B %*% b_curr)
  sigma2 <- mean(colSums((q - qstar) ^ 2)) / M

  # storage
  keep <- seq(burnin, iter, by = thin)
  nsamp <- length(keep)
  v_coef <- array(0, dim = c(ncoef, N, nsamp))
  qstar_mat <- matrix(0, M, nsamp)
  sigma_tr <- numeric(iter)
  logl_tr <- numeric(iter)
  accept <- matrix(0, N, iter)
  eps_tr <- matrix(0, N, iter)
  epsilon <- rep(sqrt(h), N)
  da <- dual_avg_init(epsilon)
  sigma_tr[1] <- sqrt(sigma2)
  eps_tr[, 1] <- epsilon
  ks <- 1

  pb <- NULL
  if (verbose)
    pb <- cli::cli_progress_bar("Running MCMC", total = iter - 1,
                                .envir = environment())

  for (m in 2:iter) {
    # update warping coefficients: the curves are conditionally independent
    # given q* and sigma^2, so all are updated in one batched HMC step
    upd <- group_hmc_update(c_curr, Phi, rep, qstar, q, sigma2, C, epsilon, L)
    c_curr <- upd$c
    accept[, m] <- upd$accept
    qtil <- upd$aux[-1, , drop = FALSE]

    # dual-averaging step-size adaptation during burn-in; keep the averaged
    # step sizes afterwards
    if (adapt && m <= burnin) {
      da_out <- dual_avg_step(da, upd$alpha, target_accept)
      epsilon <- da_out$epsilon
      da <- da_out$state
      if (m == burnin)
        epsilon <- exp(da$log_eps_bar)
    }

    # template: conjugate Gaussian draw
    P <- N * BtB / sigma2 + diag(nbasis) / qpriorvar
    R <- chol((P + t(P)) / 2)
    rhs <- as.numeric(t(B) %*% rowSums(qtil) / sigma2)
    mu <- backsolve(R, backsolve(R, rhs, transpose = TRUE))
    b_curr <- mu + backsolve(R, stats::rnorm(nbasis))
    qstar <- as.numeric(B %*% b_curr)

    # noise variance: conjugate inverse gamma draw
    SSE <- colSums((qtil - qstar) ^ 2)
    sigma2 <- 1 / stats::rgamma(1, shape = N * M / 2 + alpha0,
                                scale = 1 / (sum(SSE) / 2 + beta0))

    # centering constraint (Karcher mean of warps = identity): project onto it
    # every ncenter iterations
    if (ncenter > 0 && (m < burnin || center_all) && m %% ncenter == 0) {
      gam_c <- bayes_coef_to_gam(c_curr, Phi, rep)$gam
      gam_inv <- karcher_mean_inverse(gam_c, t01)
      for (ii in 1:N)
        gam_c[, ii] <- stats::approx(t01, gam_c[, ii], xout = gam_inv)$y
      c_curr <- bayes_gam_to_coef(gam_c, Phi, rep)
      qstar <- warp_q_gamma(qstar, t01, gam_inv)
      b_curr <- as.numeric(Bpinv %*% qstar)
      qstar <- as.numeric(B %*% b_curr)
    }

    sigma_tr[m] <- sqrt(sigma2)
    logl_tr[m] <- -sum(SSE) / (2 * sigma2) - N * M / 2 * log(2 * pi * sigma2)
    eps_tr[, m] <- epsilon
    if (ks <= nsamp && m == keep[ks]) {
      v_coef[, , ks] <- c_curr
      qstar_mat[, ks] <- qstar
      ks <- ks + 1
    }
    if (!is.null(pb)) cli::cli_progress_update(id = pb)
  }
  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  # posterior summaries of the warping functions
  gamma_mat <- array(0, dim = c(M, nsamp, N))
  gam_mean <- matrix(0, M, N)
  psi_mean <- matrix(0, M, N)
  for (ii in 1:N) {
    cg <- bayes_coef_to_gam(matrix(v_coef[, ii, ], ncoef, nsamp), Phi, rep)
    gamma_mat[, , ii] <- apply(cg$gam, 2,
                               function(g) stats::approx(t01, g, xout = tn)$y)
    pm <- rowMeans(cg$psi)
    pm <- pm / sqrt(trapz(t01, pm ^ 2))
    g <- cumtrapz(t01, pm ^ 2)
    g <- (g - g[1]) / (g[M] - g[1])
    gam_mean[, ii] <- stats::approx(t01, g, xout = tn)$y
    psi_mean[, ii] <- stats::approx(t01, pm, xout = tn)$y
  }

  fn <- matrix(0, M, N)
  for (ii in 1:N)
    fn[, ii] <- warp_f_gamma(f[, ii], time, gam_mean[, ii])
  qn <- f_to_srvf(fn, time)
  q0 <- f_to_srvf(f, time)

  # template on the original time scale (SRVF scales by 1/sqrt(T))
  mq01 <- rowMeans(qstar_mat)
  mqn <- stats::approx(t01, mq01, xout = tn)$y / sqrt(time[M] - time[1])
  fmean <- mean(f[1, ]) + spline_cumintegral(time, mqn * abs(mqn))

  std_f0 <- apply(f, 1, stats::sd)
  std_fn <- apply(fn, 1, stats::sd)
  fgam <- matrix(0, M, N)
  for (ii in 1:N)
    fgam[, ii] <- warp_f_gamma(fmean, time, gam_mean[, ii])

  mcmc <- list(
    v_coef = v_coef,
    sigma = sigma_tr,
    logl = logl_tr,
    accept = accept[, -1, drop = FALSE],
    epsilon = eps_tr,
    warp_rep = warp_rep,
    burnin = burnin,
    iter = iter
  )
  if (extrainfo) {
    mcmc$gamma_mat <- gamma_mat
    gstats <- array(0, dim = c(M, 2, N))
    for (ii in 1:N)
      gstats[, , ii] <- t(statsFun_gam(gamma_mat[, , ii]))
    mcmc$gamma_stats <- gstats
    mcmc$q_star_mat <- apply(qstar_mat, 2,
                             function(qq) stats::approx(t01, qq, xout = tn)$y) /
      sqrt(time[M] - time[1])
    mcmc$q_star_stats <- t(statsFun_gam(mcmc$q_star_mat))
    mcmc$fmean_mat <- srvf_to_f(mcmc$q_star_mat, time, rep(mean(f[1, ]), nsamp))
    mcmc$fmean_stats <- t(statsFun_gam(mcmc$fmean_mat))
  }

  out <- list(
    time = time,
    f0 = f,
    q0 = q0,
    fn = fn,
    qn = qn,
    fmean = fmean,
    mqn = mqn,
    warping_functions = gam_mean,
    psi = psi_mean,
    original_variance = trapz(time, std_f0 ^ 2),
    amplitude_variance = trapz(time, std_fn ^ 2),
    phase_variance = trapz(time, apply(fgam, 1, stats::var)),
    qun = logl_tr[-1],
    inverse_average_warping_function = invertGamma(SqrtMean(gam_mean)$gam_mu),
    rsamps = FALSE,
    mcmc = mcmc,
    call = list(
      iter = iter,
      burnin = burnin,
      thin = thin,
      ncoef = ncoef,
      basis = basis,
      warp_rep = warp_rep,
      nbasis = nbasis,
      init = init,
      ncenter = ncenter,
      center_all = center_all,
      centroid_type = "mean",
      method = "bayesian_infHMC"
    )
  )

  class(out) <- "fdawarp"

  out
}

# geometry of the group alignment likelihood for a batch of curves (one per
# column): aux = [SSE; q_i warped]
group_geom <- function(c, Phi, rep, qstar, Q, sigma2) {
  o <- bayes_warp_geom(c, Phi, rep, qstar, Q, sigma2)
  list(nll = as.numeric(o$nll), grad = o$grad,
       aux = rbind(as.numeric(o$SSE), o$q2w))
}

# one batched infinite-dimensional HMC update of the warps of curves Q
group_hmc_update <- function(c, Phi, rep, qstar, Q, sigma2, C, epsilon, L) {
  geom <- function(x) group_geom(x, Phi, rep, qstar, Q, sigma2)
  g0 <- geom(c)
  infhmc_step_batch(c, g0$nll, g0$grad, g0$aux, geom, C, epsilon, L)
}
