# Bayesian HMC alignment (pairwise and group). Mirrors BayesHMCTest.m in the
# MATLAB package: the exact likelihood gradient is checked against finite
# differences for both warp representations, and short pairwise and group runs
# must return valid warps close to the dynamic programming solution.

small_problem <- function(TT = 101L, ncoef = 10L) {
  t <- seq(0, 1, length.out = TT)
  list(
    t = t,
    Phi = fdasrvf:::bayes_warp_basis(t, ncoef, "fourier"),
    q1 = f_to_srvf(exp(-(t - 0.4) ^ 2 / 0.01), t),
    q2 = f_to_srvf(exp(-(t - 0.6) ^ 2 / 0.01), t)
  )
}

test_that("the warp likelihood gradient matches finite differences", {
  sp <- small_problem()
  withr::local_seed(1)
  cc <- 0.2 * stats::rnorm(ncol(sp$Phi))
  for (rep in 0:1) {
    g <- as.numeric(fdasrvf:::bayes_warp_geom(cc, sp$Phi, rep, sp$q1, sp$q2, 0.1)$grad)
    e <- 1e-7
    fd <- vapply(seq_along(cc), function(i) {
      cp <- cc; cp[i] <- cp[i] + e
      cm <- cc; cm[i] <- cm[i] - e
      (fdasrvf:::bayes_warp_geom(cp, sp$Phi, rep, sp$q1, sp$q2, 0.1)$nll -
         fdasrvf:::bayes_warp_geom(cm, sp$Phi, rep, sp$q1, sp$q2, 0.1)$nll) / (2 * e)
    }, numeric(1))
    expect_lt(max(abs(fd - g)) / max(abs(fd)), 1e-5)
  }
})

test_that("the shooting representation rejects warps outside the injectivity radius", {
  sp <- small_problem()
  cc <- rep(0, ncol(sp$Phi))
  cc[1] <- 4 # |v| = 4 > pi
  expect_equal(fdasrvf:::bayes_warp_geom(cc, sp$Phi, 0, sp$q1, sp$q2, 0.1)$nll[1], Inf)
})

test_that("coefficients and warping functions round trip", {
  sp <- small_problem()
  TT <- nrow(sp$Phi)
  withr::local_seed(2)
  cc <- matrix(0.1 * stats::rnorm(3 * ncol(sp$Phi)), ncol(sp$Phi), 3)
  for (rep in 0:1) {
    cg <- fdasrvf:::bayes_coef_to_gam(cc, sp$Phi, rep)
    expect_equal(cg$gam[1, ], rep(0, 3))
    expect_equal(cg$gam[TT, ], rep(1, 3))
    expect_true(all(apply(cg$gam, 2, diff) >= 0))
    expect_lt(max(abs(fdasrvf:::bayes_gam_to_coef(cg$gam, sp$Phi, rep) - cc)), 1e-2)
  }
})

test_that("the batched warp geometry matches the single-warp one", {
  sp <- small_problem()
  withr::local_seed(6)
  Q <- f_to_srvf(simu_data$f[, 1:5], sp$t)
  cc <- matrix(0.2 * stats::rnorm(5 * ncol(sp$Phi)), ncol(sp$Phi), 5)
  for (rep in 0:1) {
    o <- fdasrvf:::bayes_warp_geom(cc, sp$Phi, rep, sp$q1, Q, 0.1)
    for (i in 1:5) {
      o1 <- fdasrvf:::bayes_warp_geom(cc[, i], sp$Phi, rep, sp$q1, Q[, i], 0.1)
      expect_equal(o$nll[i], o1$nll[1])
      expect_equal(o$SSE[i], o1$SSE[1])
      expect_equal(o$grad[, i], as.numeric(o1$grad))
      expect_equal(o$gam[, i], as.numeric(o1$gam))
      expect_equal(o$q2w[, i], as.numeric(o1$q2w))
    }
  }
})

test_that("the batched HMC step samples the Gaussian prior", {
  # nll = 0: the target is the prior N(0, diag(C)), which the rotation
  # integrates exactly (every proposal is accepted)
  withr::local_seed(8)
  C <- c(1, 0.25, 0.04)
  zero_geom <- function(c)
    list(nll = rep(0, ncol(c)), grad = matrix(0, nrow(c), ncol(c)), aux = NULL)
  cc <- matrix(0, 3, 2000)
  nll <- rep(0, 2000)
  g <- matrix(0, 3, 2000)
  acc <- NULL
  for (k in 1:20) {
    o <- fdasrvf:::infhmc_step_batch(cc, nll, g, NULL, zero_geom, C, 0.7, 5)
    cc <- o$c; nll <- o$nll; g <- o$g; acc <- o$accept
  }
  expect_true(all(acc))
  expect_equal(apply(cc, 1, stats::var), C, tolerance = 0.1)
})

test_that("pairwise Bayesian HMC alignment recovers the DP warp", {
  skip_on_cran()
  withr::local_seed(4)
  f <- simu_data$f
  tt <- simu_data$time
  y1 <- f[, 1] + 0.01 * stats::rnorm(length(tt))
  y2 <- f[, 21] + 0.01 * stats::rnorm(length(tt))
  gdp <- optimum.reparam(f_to_srvf(f[, 1], tt), tt, f_to_srvf(f[, 21], tt), tt)
  for (wr in c("shooting", "clr")) {
    out <- pair_align_functions_infHMC(y1, y2, tt, iter = 1000, thin = 1,
                                       L = 10, warp_rep = wr,
                                       extrainfo = TRUE, verbose = FALSE)
    na <- length(out$theta_accept)
    expect_gt(mean(out$theta_accept[500:na]), 0.3)
    expect_lt(max(abs(out$gamma - gdp)), 0.08)
    expect_equal(dim(out$gamma_mat), c(length(tt), 501L))
  }
})

test_that("group Bayesian HMC alignment recovers the DP warps", {
  skip_on_cran()
  withr::local_seed(5)
  f <- simu_data$f
  tt <- simu_data$time
  M <- length(tt)
  N <- ncol(f)
  d <- suppressMessages(time_warping(f, tt, lambda = 0, max_iter = 10L))
  out <- time_warping_bayes_infHMC(f, tt, iter = 600, verbose = FALSE)
  expect_s3_class(out, "fdawarp")
  expect_equal(dim(out$warping_functions), c(M, N))
  expect_true(all(apply(out$warping_functions, 2, diff) >= -1e-12))
  expect_lt(mean(abs(out$warping_functions - d$warping_functions)), 0.02)
  nacc <- ncol(out$mcmc$accept)
  expect_gt(mean(out$mcmc$accept[, 300:nacc]), 0.4)
  fn <- vapply(1:N, function(k)
    warp_f_gamma(f[, k], tt, out$warping_functions[, k]), numeric(M))
  expect_lt(max(abs(fn - out$fn)), 1e-12)
  expect_equal(dim(out$mcmc$gamma_stats), c(M, 2L, N))
})

test_that("group Bayesian pCN alignment reduces the amplitude variance", {
  skip_on_cran()
  withr::local_seed(11)
  f <- simu_data$f
  tt <- simu_data$time
  out <- time_warping_bayes(f, tt, iter = 400, burnin = 200, verbose = FALSE)
  expect_s3_class(out, "fdawarp")
  expect_equal(dim(out$warping_functions), c(length(tt), ncol(f)))
  expect_true(all(apply(out$warping_functions, 2, diff) >= -1e-10))
  expect_lt(out$amplitude_variance, out$original_variance)
  expect_equal(dim(out$mcmc$gamma_mat), c(length(tt), 201L, ncol(f)))
})
