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
