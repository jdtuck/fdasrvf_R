q <- f_to_srvf(simu_data$f, simu_data$time)
time <- simu_data$time
pens <- c("none", "roughness", "l2gam", "l2psi", "geodesic")

test_that("`optimum.reparam()` returns a warping function for every penalty", {
  for (method in c("DP", "DPo", "RBFGS")) {
    for (pen in pens) {
      gam <- optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0.01,
                             pen = pen, method = method)
      expect_length(gam, length(time))
      expect_equal(gam[1], 0)
      expect_equal(gam[length(gam)], 1)
      expect_true(all(diff(gam) >= -1e-12))
    }
  }
})

test_that("`optimum.reparam()` penalty has no effect when lambda is 0", {
  for (method in c("DP", "DPo", "RBFGS")) {
    gam0 <- optimum.reparam(q[, 1], time, q[, 2], time, method = method)
    for (pen in pens) {
      gam <- optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0,
                             pen = pen, method = method)
      expect_equal(gam, gam0)
    }
  }
})

test_that("`optimum.reparam()` penalty 'none' ignores lambda", {
  for (method in c("DP", "DPo", "RBFGS")) {
    gam0 <- optimum.reparam(q[, 1], time, q[, 2], time, method = method)
    gam <- optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0.5,
                           pen = "none", method = method)
    expect_equal(gam, gam0)
  }
})

test_that("`optimum.reparam()` passes the penalty to the DP method", {
  gam <- lapply(c("none", "roughness", "l2gam"), function(pen) {
    optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0.1, pen = pen,
                    method = "DP")
  })
  expect_false(isTRUE(all.equal(gam[[1]], gam[[2]])))
  expect_false(isTRUE(all.equal(gam[[2]], gam[[3]])))
})

test_that("`optimum.reparam()` rejects an unknown penalty", {
  expect_error(
    optimum.reparam(q[, 1], time, q[, 2], time, pen = "bogus"),
    "invalid penalty selection"
  )
})

# f2 = f1 o gam_true, so aligning q2 to q1 should recover the inverse of
# gam_true (as in MATLAB tests/DynamicProgrammingTest.m)
kw_time <- seq(0, 1, length.out = 101)
kw_f <- function(s) sin(3 * pi * s) + 0.5 * cos(7 * pi * s) + 2 * s
kw_gam_true <- kw_time + 0.5 * sin(pi * kw_time) / pi
kw_q1 <- f_to_srvf(kw_f(kw_time), kw_time)
kw_q2 <- f_to_srvf(kw_f(kw_gam_true), kw_time)
kw_gam_inv <- stats::approx(kw_gam_true, kw_time, xout = kw_time)$y
kw_error <- function(gam) {
  qw <- warp_q_gamma(kw_q2, kw_time, gam)
  sqrt(trapz(kw_time, (kw_q1 - qw)^2))
}

test_that("`optimum.reparam()` DPo recovers a known warp", {
  # rounding q2 to the nearest grid point at scl = 1 left the warp up to
  # 0.63 grid cells off
  gam <- optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = "DPo")
  expect_lt(max(abs(gam - kw_gam_inv)) * (length(kw_time) - 1), 0.5)
})

test_that("`optimum.reparam()` DP and DPo fit as well as the true warp", {
  # the old DPo cost fit 1.35x and the old DP cost (piecewise-constant SRVFs,
  # nbhd_dim = 7) 1.52x worse than the true warp
  for (method in c("DP", "DPo")) {
    gam <- optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = method)
    expect_lt(kw_error(gam), 1.1 * kw_error(kw_gam_inv))
  }
})

test_that("`optimum.reparam()` passes scl to the DPo solver", {
  gam <- optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = "DPo")
  expect_identical(
    optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = "DPo", scl = 3),
    gam
  )
  gam1 <- optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = "DPo",
                          scl = 1)
  expect_length(gam1, length(gam))
  expect_true(all(diff(gam1) >= 0))
  expect_false(isTRUE(all.equal(gam1, gam)))
  expect_error(
    optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = "DPo", scl = 0),
    "scl"
  )
})

test_that("`optimum.reparam()` RBFGS fits as well as the true warp", {
  # "RBFGS" refines a coarse "DP" warp on SRVFs resampled scl times finer;
  # started from the identity it stopped half a grid cell off
  gam <- optimum.reparam(kw_q1, kw_time, kw_q2, kw_time, method = "RBFGS")
  expect_lt(max(abs(gam - kw_gam_inv)) * (length(kw_time) - 1), 0.5)
  expect_lt(kw_error(gam), 1.1 * kw_error(kw_gam_inv))
})

test_that("`optimum.reparam()` RBFGS distances match DP", {
  # on random pairs the identity start stopped in local minima, up to 3x the
  # "DP" distance
  bumps <- function(t) {
    f <- rep(0, length(t))
    for (j in 1:3)
      f <- f + stats::rnorm(1) *
        exp(-(t - stats::runif(1))^2 / (2 * (0.05 + 0.1 * stats::runif(1))^2))
    f
  }
  amp_dist <- function(f1, f2, method) {
    q1 <- f_to_srvf(f1, kw_time)
    q2 <- f_to_srvf(f2, kw_time)
    gam <- optimum.reparam(q1, kw_time, q2, kw_time, method = method)
    qw <- warp_q_gamma(q2, kw_time, gam)
    sqrt(trapz(kw_time, (q1 - qw)^2))
  }
  withr::local_seed(0)
  r <- replicate(10, {
    f1 <- bumps(kw_time)
    f2 <- bumps(kw_time)
    amp_dist(f1, f2, "RBFGS") / amp_dist(f1, f2, "DP")
  })
  expect_lt(stats::median(r), 1.0)
  # MATLAB bounds the max by 1.1 on its random pairs; on R's, a warp steeper
  # than the 101 point grid resolves ties DP on the scl finer grid it is
  # fitted on but measures up to 1.2x here
  expect_lt(max(r), 1.25)
})

test_that("`rlbfgs()` starts from gam0 and reports its stop", {
  q1n <- kw_q1 / sqrt(sum(kw_q1^2))
  q2n <- kw_q2 / sqrt(sum(kw_q2^2))
  out <- rlbfgs(q1n, q2n, kw_time, 30, 0, 0, kw_gam_inv)
  expect_length(out$info, 4)
  expect_true(out$info[2] %in% 1:3)
  expect_true(all(diff(out$gam) >= 0))
  expect_lt(max(abs(out$gam - kw_gam_inv)) * (length(kw_time) - 1), 0.5)
  expect_error(rlbfgs(q1n, q2n, kw_time, 30, 0, 0, kw_gam_inv[-1]), "gam0")
})
