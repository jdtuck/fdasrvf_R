# Tests for the curve/function <-> SRVF representation converters:
#   * discrete2curve / curve2srvf / srvf2curve (functional form, R/refactoring.R),
#   * f_to_srvf / srvf_to_f (discrete uni-dimensional form),
#   * warp_curve / warp_srvf under the identity warping.
# Correctness is anchored on round-trip identities. Note that SRVF/curve
# round-trips recover the curve only up to an additive constant (the SRVF
# discards absolute position), so curves are compared after centering.

test_that("`srvf2curve(curve2srvf())` recovers a curve up to translation", {
  beta1 <- fdasrvf::beta[, , 1, 1]
  q <- curve2srvf(beta1)
  betafun_rec <- srvf2curve(q)
  s <- seq(0, 1, length.out = ncol(beta1))
  beta_rec <- betafun_rec(s)          # L x M matrix
  # Compare centered curves (SRVF is invariant to translation).
  center <- function(x) x - rowMeans(x)
  expect_equal(dim(beta_rec), dim(beta1))
  expect_equal(center(beta_rec), center(beta1), tolerance = 1e-2)
})

test_that("`srvf_to_f(f_to_srvf())` round-trips a 1-D function", {
  f <- simu_data$f[, 1]
  t <- simu_data$time
  q <- f_to_srvf(f, t)
  expect_length(q, length(f))
  # `cumtrapz()` inside srvf_to_f() returns an M x 1 matrix; drop dims.
  f_rec <- as.numeric(srvf_to_f(q, t, f0 = f[1]))
  expect_length(f_rec, length(f))
  expect_equal(f_rec, f, tolerance = 1e-2)
})

test_that("`warp_curve()` under the identity warping is a no-op", {
  betafun <- discrete2curve(fdasrvf::beta[, , 1, 1])
  warped <- warp_curve(betafun, get_identity_warping())
  s <- seq(0, 1, length.out = 25)
  expect_equal(warped(s), betafun(s), tolerance = 1e-6)
})

test_that("`warp_srvf()` under the identity warping is a no-op", {
  q <- curve2srvf(fdasrvf::beta[, , 1, 1])
  warped <- warp_srvf(q, get_identity_warping())
  s <- seq(0, 1, length.out = 25)
  expect_equal(warped(s), q(s), tolerance = 1e-6)
})
