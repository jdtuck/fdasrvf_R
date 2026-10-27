# Tests for the warping-function <-> tangent-space / Hilbert-sphere conversion
# helpers defined in R/geometry.R and the gamma inversion helper in
# R/invertGamma.R. Correctness is anchored on mathematical invariants:
#   * round-trip identities (converting there and back recovers the input),
#   * normalization of warping functions (start at 0, end at 1, monotone),
#   * behavior at the identity warping.
#
# Discrete gradient / cumulative-trapezoid operations introduce O(1/M) error,
# so round-trip tolerances are intentionally loose.

is_valid_warping <- function(g, tol = 1e-6) {
  # A warping function on [0, 1] must start at 0, end at 1 and be
  # monotonically non-decreasing.
  g[1] >= -tol &&
    abs(g[length(g)] - 1) < 1e-3 &&
    all(diff(g) >= -tol)
}

test_that("`gam_to_v()` / `v_to_gam()` round-trip and identity work", {
  gam <- toy_warp$gam[, 1]
  v <- gam_to_v(gam, smooth = FALSE)
  # `cumtrapz()` returns an M x 1 matrix, so drop dims before comparing values.
  gam_rec <- as.numeric(v_to_gam(v))
  expect_length(v, length(gam))
  expect_length(gam_rec, length(gam))
  expect_true(is_valid_warping(gam_rec))
  expect_equal(gam_rec, gam, tolerance = 1e-2)

  # The identity warping maps to the zero shooting vector.
  gamid <- seq(0, 1, length.out = 101)
  vid <- gam_to_v(gamid, smooth = FALSE)
  expect_equal(vid, rep(0, 101), tolerance = 1e-2)
})

test_that("`gam_to_psi()` / `psi_to_gam()` round-trip and identity work", {
  gam <- toy_warp$gam[, 1]
  psi <- gam_to_psi(gam, smooth = FALSE)
  gam_rec <- as.numeric(psi_to_gam(psi))
  expect_length(psi, length(gam))
  expect_true(is_valid_warping(gam_rec))
  expect_equal(gam_rec, gam, tolerance = 1e-2)

  # The identity warping maps to the constant-one function on the sphere.
  gamid <- seq(0, 1, length.out = 101)
  psiid <- gam_to_psi(gamid, smooth = FALSE)
  expect_equal(psiid, rep(1, 101), tolerance = 1e-2)
})

test_that("`gam_to_h()` / `h_to_gam()` round-trip works", {
  gam <- toy_warp$gam[, 1]
  h <- gam_to_h(gam, smooth = FALSE)
  gam_rec <- as.numeric(h_to_gam(h))
  expect_length(h, length(gam))
  expect_true(is_valid_warping(gam_rec))
  expect_equal(gam_rec, gam, tolerance = 1e-2)
})

test_that("conversion helpers preserve the shape of matrix input", {
  gam <- toy_warp$gam
  d <- dim(gam)
  expect_equal(dim(gam_to_v(gam, smooth = FALSE)), d)
  expect_equal(dim(gam_to_psi(gam, smooth = FALSE)), d)
  expect_equal(dim(gam_to_h(gam, smooth = FALSE)), d)
  expect_equal(dim(v_to_gam(gam_to_v(gam, smooth = FALSE))), d)
  expect_equal(dim(psi_to_gam(gam_to_psi(gam, smooth = FALSE))), d)
  expect_equal(dim(h_to_gam(gam_to_h(gam, smooth = FALSE))), d)
})

test_that("`invertGamma()` inverts a warping function", {
  gam <- toy_warp$gam[, 1]
  gamI <- invertGamma(gam)
  expect_length(gamI, length(gam))
  expect_true(is_valid_warping(gamI))
  # Inverting twice recovers the original warping function.
  expect_equal(invertGamma(gamI), gam, tolerance = 1e-2)

  # The identity warping is its own inverse.
  gamid <- seq(0, 1, length.out = 101)
  expect_equal(invertGamma(gamid), gamid, tolerance = 1e-2)
})

test_that("`inv_exp_map()` is accurate for nearby points", {
  # acos of the inner product would round the angle to zero or ~1e-8
  M <- 101
  s <- seq(0, 1, length.out = M)
  psi <- rep(1, M)
  v0 <- sin(2 * pi * s)
  v <- v0 - trapz(s, v0 * psi) * psi
  v <- 1e-7 * v / sqrt(trapz(s, v^2))
  p2 <- cos(1e-7) * psi + sin(1e-7) * v / 1e-7
  v2 <- inv_exp_map(psi, p2)
  expect_equal(sqrt(trapz(s, v2^2)), 1e-7, tolerance = 1e-6)
})

test_that("`gam_to_psi()` lies on the Hilbert sphere", {
  s <- seq(0, 1, length.out = 101)
  gam <- cbind(s, (exp(3 * s) - 1) / (exp(3) - 1))
  psi <- gam_to_psi(gam)
  expect_equal(apply(psi, 2, function(p) sqrt(trapz(s, p^2))), c(1, 1),
               tolerance = 1e-12)
})

test_that("`gam_to_h()` is finite for flat warp segments", {
  s <- seq(0, 1, length.out = 101)
  gam <- stats::approx(c(0, 0.3, 0.5, 1), c(0, 0.4, 0.4, 1), xout = s)$y
  expect_true(all(is.finite(gam_to_h(gam))))
})

test_that("`elastic.distance()` of a function to itself is zero", {
  t <- seq(0, 1, length.out = 101)
  f <- sin(2 * pi * t) + t
  out <- elastic.distance(f, f, t)
  expect_lt(out$Dx, 1e-10)
})
