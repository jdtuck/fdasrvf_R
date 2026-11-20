# Group-wise Bayesian alignment with the Z-mixture pCN algorithm.
#
# Ported from the MATLAB package (fdawarp.time_warping_bayes).

# B-spline basis in the list(x, matrix) form used by the expomap helpers
basis.bspline <- function(f.domain, numBasis, df = 4) {
  list(x = f.domain, matrix = create_basismatrix(f.domain, numBasis, df))
}

# SSE of the group alignment: for each curve, q[, ii] is warped by g[, ii]
# and compared with the template q_star
f.SSEg.group <- function(g, q, q_star) {
  obs.domain <- q$x
  N <- ncol(q$y)
  out <- numeric(N)
  for (ii in 1:N) {
    g_tmp <- list(x = g$x, y = g$y[, ii])
    exp1g.temp <- f.predictfunction(f = f.exp1(g_tmp), at = obs.domain)
    pt <- c(0, cuL2norm2(x = obs.domain, y = exp1g.temp$y))
    q_tmp <- list(x = obs.domain, y = q$y[, ii])
    tmp <- f.predictfunction(f = q_tmp, at = pt)
    out[ii] <- sum((tmp$y * exp1g.temp$y - q_star) ^ 2)
  }
  out
}

f.logl.group <- function(g, q, q_star, var1, SSEg = 0) {
  if (all(SSEg == 0))
    SSEg <- f.SSEg.group(g = g, q = q, q_star = q_star)
  n <- length(q$x)
  sum(n * log(1 / sqrt(2 * pi)) - n * log(sqrt(var1)) - SSEg / (2 * var1))
}

# one pCN update of the warping coefficients of a single curve
f.updateg.group <- function(g.coef.curr, g.basis, var1.curr, q, q_star,
                            propose.g.coef) {
  g.coef.prop <- propose.g.coef(g.coef.curr = g.coef.curr)
  tst <- f.exp1(f.basistofunction(f.domain = g.basis$x,
                                  coef = g.coef.prop$prop, basis = g.basis))
  while (min(tst$y) < 0) {
    g.coef.prop <- propose.g.coef(g.coef.curr = g.coef.curr)
    tst <- f.exp1(f.basistofunction(f.domain = g.basis$x,
                                    coef = g.coef.prop$prop, basis = g.basis))
  }

  gc <- list(x = g.basis$x,
             y = matrix(f.basistofunction(f.domain = g.basis$x,
                                          coef = g.coef.curr,
                                          basis = g.basis)$y, ncol = 1))
  gp <- list(x = g.basis$x,
             y = matrix(f.basistofunction(f.domain = g.basis$x,
                                          coef = g.coef.prop$prop,
                                          basis = g.basis)$y, ncol = 1))
  SSE.curr <- f.SSEg.group(g = gc, q = q, q_star = q_star)
  SSE.prop <- f.SSEg.group(g = gp, q = q, q_star = q_star)

  ratio <- min(1, exp(-(SSE.prop - SSE.curr) / (2 * var1.curr)))
  if (stats::runif(1) <= ratio)
    list(g.coef = g.coef.prop$prop, accept = TRUE, zpcnInd = g.coef.prop$ind)
  else
    list(g.coef = g.coef.curr, accept = FALSE, zpcnInd = g.coef.prop$ind)
}

# one pCN update of the template coefficients
f.updateq.group <- function(g.coef.curr, g.basis, var1.curr, q,
                            q.star.coef.curr, q.basis, SSE.curr,
                            propose.q.star) {
  q.star.coef.prop <- propose.q.star(q.star.coef.curr)
  q.star.prop <- f.basistofunction(f.domain = q.basis$x,
                                   coef = q.star.coef.prop$prop,
                                   basis = q.basis)
  q.star.curr <- f.basistofunction(f.domain = q.basis$x,
                                   coef = q.star.coef.curr,
                                   basis = q.basis)

  g <- list(x = g.basis$x,
            y = g.basis$matrix %*% g.coef.curr)

  if (all(SSE.curr == 0))
    SSE.curr <- f.SSEg.group(g = g, q = q, q_star = q.star.curr$y)
  SSE.prop <- f.SSEg.group(g = g, q = q, q_star = q.star.prop$y)

  logl.curr <- f.logl.group(g, q, q.star.curr$y, var1.curr, SSE.curr)
  logl.prop <- f.logl.group(g, q, q.star.prop$y, var1.curr, SSE.prop)

  ratio <- min(1, exp(logl.prop - logl.curr))
  if (stats::runif(1) <= ratio)
    list(q.star.coef = q.star.coef.prop$prop, accept = TRUE,
         zpcnInd = q.star.coef.prop$ind)
  else
    list(q.star.coef = q.star.coef.curr, accept = FALSE,
         zpcnInd = q.star.coef.prop$ind)
}

#' Bayesian group-wise alignment using the Z-mixture pCN algorithm
#'
#' This function aligns a collection of functions using the elastic
#' square-root slope (SRSF) framework with a Bayesian approach. It is based on
#' mapping warping functions to a hypersphere, and a subsequent exponential
#' mapping to a tangent space. In the tangent space, the Z-mixture pCN
#' algorithm is used to explore both local and global structure in the
#' posterior distribution.
#'
#' The Z-mixture pCN algorithm uses a mixture distribution for the proposal
#' distribution, controlled by input parameter `zpcn`. The `zpcn$betas` must
#' be between 0 and 1, and are the coefficients of the mixture components,
#' with larger coefficients corresponding to larger shifts in parameter
#' space. The `zpcn$probs` give the probability of each shift size.
#'
#' @param f matrix (\eqn{M} x \eqn{N}) of \eqn{N} functions with \eqn{M}
#'   samples
#' @param time vector of length \eqn{M} describing the sample points
#' @param iter number of iterations
#' @param burnin number of burn-in MCMC iterations; defaults to
#'   `min(5e3, floor(iter / 2))`
#' @param alpha0,beta0 inverse gamma prior parameters for the model variance
#' @param zpcn list of mixture coefficients and prior probabilities for the
#'   Z-mixture pCN algorithm, of the form `list(betas, probs)`, where `betas`
#'   and `probs` are numeric vectors of equal length
#' @param propvar variance of the proposal distribution
#' @param init.coef initial coefficients of the warping functions, a matrix
#'   with an even number of rows and \eqn{N} columns
#' @param nbasis number of B-spline basis functions used for the template
#' @param npoints number of sample interpolation points
#' @param ncenter number of iterations between centering of the warping
#'   functions (`0` to never recenter)
#' @param extrainfo whether additional information is returned
#' @param verbose whether to show a progress bar
#'
#' @return An object of class `fdawarp` which is a list with the same
#'   components as [time_warping()], where `warping_functions` holds the
#'   posterior mean warping function of each curve, together with a `mcmc`
#'   component containing
#' \item{g_coef}{array (`iter` x \eqn{D} x \eqn{N}) of posterior draws of the
#'   warping coefficients}
#' \item{sigma1}{posterior draws of the model standard deviation}
#' \item{logl}{log-likelihood trace}
#' \item{burnin, iter}{number of burn-in and total iterations}
#' \item{extrainfo}{whether the extra information was stored}
#' \item{accept}{acceptance of the template samples (if `extrainfo`)}
#' \item{betas_ind}{index of the `zpcn` mixture component used at each
#'   iteration (if `extrainfo`)}
#' \item{zpcn}{the `zpcn` mixture used (if `extrainfo`)}
#' \item{gamma_mat}{array (\eqn{M} x `nsamp` x \eqn{N}) of posterior gammas
#'   (if `extrainfo`)}
#' \item{gamma_stats}{pointwise 95% credible bands of the gammas (if
#'   `extrainfo`)}
#' \item{q_star_mat}{posterior draws of the template (if `extrainfo`)}
#' \item{q_star_stats}{pointwise 95% bands of the template (if `extrainfo`)}
#' \item{fmean_mat}{posterior draws of the mean function (if `extrainfo`)}
#' \item{fmean_stats}{pointwise 95% bands of the mean function (if
#'   `extrainfo`)}
#'
#' @keywords bayesian
#' @concept srsf alignment
#'
#' @references Lu, Y., Herbei, R., and Kurtek, S. (2017). Bayesian
#'   registration of functions with a Gaussian process prior. Journal of
#'   Computational and Graphical Statistics, 26(4), 894-904.
#'
#' @export
#' @examples
#' \dontrun{
#'   # This is an MCMC algorithm and takes a long time to run
#'   out <- time_warping_bayes(simu_data$f, simu_data$time, iter = 500)
#' }
time_warping_bayes <- function(f,
                               time,
                               iter = 2000L,
                               burnin = NULL,
                               alpha0 = 0.1,
                               beta0 = 0.1,
                               zpcn = list(betas = c(0.5, 0.5, 0.005, 0.0001),
                                           probs = c(0.1, 0.1, 0.7, 0.1)),
                               propvar = 4,
                               init.coef = NULL,
                               nbasis = 30L,
                               npoints = 200L,
                               ncenter = 100L,
                               extrainfo = TRUE,
                               verbose = TRUE) {
  f <- as.matrix(f)
  time <- as.numeric(time)
  M <- nrow(f)
  N <- ncol(f)
  if (length(time) != M)
    stop("Length of time must equal nrow(f)")
  if (length(zpcn$betas) != length(zpcn$probs))
    stop("In zpcn, betas must equal length of probs")
  if (is.null(init.coef))
    init.coef <- matrix(0, 20, N)
  init.coef <- as.matrix(init.coef)
  if (nrow(init.coef) %% 2 != 0)
    stop("Number of rows of init.coef must be even")
  if (ncol(init.coef) != N)
    stop("init.coef must have one column per function")
  if (nbasis %% 2 != 0)
    stop("nbasis must be even")
  if (is.null(burnin))
    burnin <- min(5e3, floor(iter / 2))
  if (burnin < 1 || burnin >= iter)
    stop("burnin must be between 1 and iter - 1")

  # Number of sig figs to report in gamma_mat
  SIG_GAM <- 13

  # normalize time to [0,1]
  tn <- (time - min(time)) / (max(time) - min(time))

  # parameter settings
  valid_index <- burnin:iter
  Mg <- nrow(init.coef) / 2
  domain.par <- seq(0, 1, length.out = npoints)
  g.basis <- basis.fourier(f.domain = domain.par, numBasis = Mg,
                           fourier.p = 1)
  Mq <- nbasis / 2
  q.basis <- basis.bspline(domain.par, 2 * Mq, 3)
  sigma1.ini <- 1

  # srsf transformation on the simulation grid
  f0 <- apply(f, 2, function(y) stats::approx(tn, y, xout = domain.par)$y)
  qo <- f_to_srvf(f0, domain.par)

  sigma.g <- propvar
  sigma.q <- (stats::quantile(qo, 0.95) - stats::quantile(qo, 0.05)) * 5

  propose.g.coef <- function(g.coef.curr) {
    probm <- c(0, cumsum(zpcn$probs))
    z <- stats::runif(1)
    for (i in seq_along(zpcn$betas)) {
      if (z <= probm[i + 1] && z > probm[i]) {
        g.coef.new <- stats::rnorm(Mg * 2, mean = 0,
                                   sd = sigma.g / rep(1:Mg, each = 2))
        return(list(prop = sqrt(1 - zpcn$betas[i] ^ 2) * g.coef.curr +
                      zpcn$betas[i] * g.coef.new,
                    ind = i))
      }
    }
  }

  propose.q.star <- function(q.star.coef.curr) {
    probm <- c(0, cumsum(zpcn$probs))
    sdvec <- sigma.q / rep(1:Mq, each = 2)
    z <- stats::runif(1)
    for (i in seq_along(zpcn$betas)) {
      if (z <= probm[i + 1] && z > probm[i]) {
        q.star.new <- stats::rnorm(Mq * 2, mean = 0, sd = sdvec)
        return(list(prop = sqrt(1 - zpcn$betas[i] ^ 2) * q.star.coef.curr +
                      zpcn$betas[i] * q.star.new,
                    ind = i))
      }
    }
  }

  for (ii in 1:N) {
    tmp <- f.exp1(f.basistofunction(f.domain = g.basis$x,
                                    coef = init.coef[, ii], basis = g.basis))
    if (min(tmp$y) < 0)
      stop(sprintf("Invalid initial value of g for sample %d", ii))
  }

  # result objects
  res.g.coef <- array(0, dim = c(iter, nrow(init.coef), N))
  res.q.star <- matrix(0, iter, length(g.basis$x))
  res.sigma1 <- numeric(iter)
  res.logl <- numeric(iter)
  res.SSE <- matrix(0, iter, N)
  res.accept <- numeric(iter)
  res.accept.betas <- numeric(iter)

  # init
  g.coef.curr <- init.coef
  sigma1.curr <- sigma1.ini
  dqq <- sqrt(colSums((qo - rowMeans(qo)) ^ 2))
  q.star.ini <- qo[, which.min(dqq)]
  q <- list(x = g.basis$x, y = qo)

  g.curr <- list(x = g.basis$x, y = g.basis$matrix %*% g.coef.curr)
  SSE.curr <- f.SSEg.group(g.curr, q, q.star.ini)
  logl.curr <- f.logl.group(g.curr, q, q.star.ini, sigma1.ini ^ 2, SSE.curr)

  res.g.coef[1, , ] <- init.coef
  res.sigma1[1] <- sigma1.ini
  res.q.star[1, ] <- q.star.ini
  res.SSE[1, ] <- SSE.curr
  res.logl[1] <- logl.curr

  q.star.coef.curr <- as.numeric(qr.solve(q.basis$matrix, q.star.ini))
  q.star.curr <- list(x = g.basis$x, y = q.star.ini)
  gamI_o <- seq(0, 1, length.out = length(g.basis$x))

  pb <- NULL
  if (verbose)
    pb <- cli::cli_progress_bar("Running MCMC", total = iter - 1,
                                .envir = environment())

  for (m in 2:iter) {
    # update g
    for (ii in 1:N) {
      q_tmp <- list(x = q$x, y = matrix(q$y[, ii], ncol = 1))
      upd <- f.updateg.group(g.coef.curr[, ii], g.basis, sigma1.curr ^ 2,
                             q_tmp, q.star.curr$y, propose.g.coef)
      g.coef.curr[, ii] <- upd$g.coef
    }
    g.curr <- list(x = g.basis$x, y = g.basis$matrix %*% g.coef.curr)
    SSE.curr <- f.SSEg.group(g.curr, q, q.star.curr$y)

    # update q_star
    upd <- f.updateq.group(g.coef.curr, g.basis, sigma1.curr ^ 2, q,
                           q.star.coef.curr, q.basis, SSE.curr, propose.q.star)
    q.star.coef.curr <- upd$q.star.coef
    accept <- upd$accept
    zpcnInd <- upd$zpcnInd
    q.star.curr <- f.basistofunction(f.domain = q.basis$x,
                                     coef = q.star.coef.curr,
                                     basis = q.basis)

    # update sigma1
    SSE.curr <- f.SSEg.group(g.curr, q, q.star.curr$y)
    newshape <- N * length(time) / 2 + alpha0
    newscale <- sum(SSE.curr) / 2 + beta0
    sigma1.curr <- sqrt(1 / stats::rgamma(1, shape = newshape,
                                          scale = 1 / newscale))
    logl.curr <- f.logl.group(g.curr, q, q.star.curr$y, sigma1.curr ^ 2,
                              SSE.curr)

    # center
    if (ncenter != 0 && m %% ncenter == 0) {
      post.gamma <- matrix(0, length(g.basis$x), N)
      for (ii in 1:N) {
        g.temp <- f.basistofunction(f.domain = domain.par,
                                    coef = g.coef.curr[, ii],
                                    basis = g.basis)
        post.gamma[, ii] <- norm_gam(f.phiinv(f.exp1(g.temp))$y)
      }
      gamI_o <- invertGamma(SqrtMean(post.gamma)$gam_mu)
      q.star.curr$y <- warp_q_gamma(q.star.curr$y, q.star.curr$x, gamI_o)
      for (k in 1:N) {
        post.gamma[, k] <- stats::approx(
          q.star.curr$x, post.gamma[, k],
          xout = (q.star.curr$x[length(q.star.curr$x)] - q.star.curr$x[1]) *
            gamI_o + q.star.curr$x[1])$y
        gamma_k <- list(x = g.basis$x, y = post.gamma[, k])
        yy <- f.exp1inv(f.phi(gamma_k))$y
        g.coef.curr[, k] <- qr.solve(g.basis$matrix, yy)
      }
      g.curr <- list(x = g.basis$x, y = g.basis$matrix %*% g.coef.curr)
      SSE.curr <- f.SSEg.group(g.curr, q, q.star.curr$y)
      logl.curr <- f.logl.group(g.curr, q, q.star.curr$y, sigma1.curr ^ 2,
                                SSE.curr)
    }

    # save update to results
    res.g.coef[m, , ] <- g.coef.curr
    res.q.star[m, ] <- q.star.curr$y
    res.sigma1[m] <- sigma1.curr
    res.SSE[m, ] <- SSE.curr
    if (extrainfo) {
      res.logl[m] <- logl.curr
      res.accept[m] <- accept
      res.accept.betas[m] <- zpcnInd
    }
    if (!is.null(pb)) cli::cli_progress_update(id = pb)
  }
  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  # calculate posterior mean of psi
  nv <- length(valid_index)
  psi.mat <- array(0, dim = c(length(domain.par), nv, N))
  for (k in 1:nv) {
    for (ii in 1:N) {
      g.temp <- f.basistofunction(f.domain = domain.par,
                                  coef = res.g.coef[valid_index[k], , ii],
                                  basis = g.basis)
      psi.mat[, k, ii] <- f.exp1(g.temp)$y
    }
  }

  post.psi <- matrix(0, M, N)
  post.gamma <- matrix(0, M, N)
  f_warped <- matrix(0, M, N)
  for (ii in 1:N) {
    pm <- f.psimean(domain.par, psi.mat[, , ii])
    result_i <- stats::approx(pm$x, pm$y, xout = tn, rule = 2)$y
    post.psi[, ii] <- result_i
    post.gamma[, ii] <- norm_gam(f.phiinv(list(x = tn, y = result_i))$y)
    f_warped[, ii] <- warp_f_gamma(f[, ii], time, post.gamma[, ii])
  }

  qn <- f_to_srvf(f_warped, time)
  q0 <- f_to_srvf(f, time)
  mqn_temp <- colMeans(res.q.star[valid_index, , drop = FALSE])
  mqn <- stats::approx(g.basis$x, mqn_temp, xout = tn, rule = 2)$y
  fmean <- mean(f[1, ]) + spline_cumintegral(time, mqn * abs(mqn))

  std_f0 <- apply(f, 1, stats::sd)
  std_fn <- apply(f_warped, 1, stats::sd)
  fgam <- matrix(0, M, N)
  for (ii in 1:N)
    fgam[, ii] <- warp_f_gamma(fmean, time, post.gamma[, ii])

  mcmc <- list(g_coef = res.g.coef, sigma1 = res.sigma1, logl = res.logl,
               burnin = burnin, iter = iter, extrainfo = extrainfo)
  if (extrainfo) {
    gamma_mat <- array(0, dim = c(M, nv, N))
    gamma_stats <- array(0, dim = c(M, 2, N))
    for (jj in 1:N) {
      for (ii in 1:nv) {
        result_i <- stats::approx(domain.par, psi.mat[, ii, jj], xout = tn,
                                  rule = 2)$y
        gamma_mat[, ii, jj] <- round(
          norm_gam(f.phiinv(list(x = tn, y = result_i))$y), SIG_GAM)
      }
      gamma_stats[, , jj] <- t(statsFun_gam(gamma_mat[, , jj]))
    }
    mcmc$accept <- res.accept[-1]
    mcmc$betas_ind <- res.accept.betas[-1]
    mcmc$zpcn <- zpcn
    mcmc$gamma_mat <- gamma_mat
    mcmc$gamma_stats <- gamma_stats
    mcmc$q_star_mat <- t(res.q.star[valid_index, , drop = FALSE])
    mcmc$q_star_stats <- t(statsFun_gam(mcmc$q_star_mat))
    mcmc$fmean_mat <- srvf_to_f(mcmc$q_star_mat, q.basis$x,
                                rep(mean(f[1, ]), nv))
    mcmc$fmean_stats <- t(statsFun_gam(mcmc$fmean_mat))
  }

  out <- list(
    time = time,
    f0 = f,
    q0 = q0,
    fn = f_warped,
    qn = qn,
    fmean = fmean,
    mqn = mqn,
    warping_functions = post.gamma,
    psi = post.psi,
    original_variance = trapz(time, std_f0 ^ 2),
    amplitude_variance = trapz(time, std_fn ^ 2),
    phase_variance = trapz(time, apply(fgam, 1, stats::var)),
    qun = res.logl,
    inverse_average_warping_function =
      stats::approx(g.basis$x, gamI_o, xout = tn, rule = 2)$y,
    rsamps = FALSE,
    mcmc = mcmc,
    call = list(
      iter = iter,
      burnin = burnin,
      nbasis = nbasis,
      npoints = npoints,
      ncenter = ncenter,
      centroid_type = "mean",
      method = "bayesian"
    )
  )

  class(out) <- "fdawarp"

  out
}
