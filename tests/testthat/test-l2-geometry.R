# Tests for the L^2 and Hilbert-sphere geometry helpers defined in
# R/refactoring.R. SRVFs are built from the sample curves in `beta` using the
# documented `curve2srvf()` / `discrete2curve()` constructors. Correctness is
# anchored on the defining properties of an inner-product space and a metric:
#   * ||q|| = sqrt(<q, q>),
#   * distance from a point to itself is 0,
#   * symmetry of distance and inner product,
#   * projection onto the sphere yields unit norm,
#   * geodesic distances live in [0, pi].

test_that("`get_l2_norm()` and `get_l2_inner_product()` are consistent", {
  q <- curve2srvf(fdasrvf::beta[, , 1, 1])
  nrm <- get_l2_norm(q)
  ip <- get_l2_inner_product(q, q)
  expect_true(nrm > 0)
  expect_equal(nrm, sqrt(ip), tolerance = 1e-6)
})

test_that("`get_l2_inner_product()` is symmetric", {
  q1 <- curve2srvf(fdasrvf::beta[, , 1, 1])
  q2 <- curve2srvf(fdasrvf::beta[, , 1, 2])
  expect_equal(
    get_l2_inner_product(q1, q2),
    get_l2_inner_product(q2, q1),
    tolerance = 1e-6
  )
})

test_that("`get_l2_distance()` is a metric (identity + symmetry)", {
  q1 <- curve2srvf(fdasrvf::beta[, , 1, 1])
  q2 <- curve2srvf(fdasrvf::beta[, , 1, 2])
  expect_equal(get_l2_distance(q1, q1), 0, tolerance = 1e-6)
  d12 <- get_l2_distance(q1, q2)
  d21 <- get_l2_distance(q2, q1)
  expect_true(d12 > 0)
  expect_equal(d12, d21, tolerance = 1e-6)
})

test_that("`to_hilbert_sphere()` produces a unit-norm SRVF", {
  q <- curve2srvf(fdasrvf::beta[, , 1, 1])
  qh <- to_hilbert_sphere(q)
  expect_equal(get_l2_norm(qh), 1, tolerance = 1e-6)
})

test_that("`get_hilbert_sphere_distance()` is a valid geodesic distance", {
  q1 <- curve2srvf(fdasrvf::beta[, , 1, 1])
  q2 <- curve2srvf(fdasrvf::beta[, , 1, 2])
  q1h <- to_hilbert_sphere(q1)
  q2h <- to_hilbert_sphere(q2)
  expect_equal(get_hilbert_sphere_distance(q1h, q1h), 0, tolerance = 1e-6)
  d <- get_hilbert_sphere_distance(q1h, q2h)
  expect_true(d >= 0 && d <= pi)
})

test_that("`get_identity_warping()` returns the identity and its derivatives", {
  id <- get_identity_warping()
  s <- seq(0, 1, length.out = 11)
  expect_equal(id(s), s)
  expect_equal(id(s, deriv = 1), rep(1, length(s)))
  expect_equal(id(s, deriv = 2), rep(0, length(s)))
})

test_that("`get_warping_distance()` is zero to itself and non-negative", {
  gam1 <- discrete2warping(toy_warp$gam[, 1])
  gam2 <- discrete2warping(toy_warp$gam[, 2])
  expect_equal(get_warping_distance(gam1, gam1), 0, tolerance = 1e-6)
  d <- get_warping_distance(gam1, gam2)
  expect_true(d >= 0)
  expect_true(get_warping_distance(gam1, get_identity_warping()) >= 0)
})

test_that("`get_curve_centroid()` returns a length-L vector", {
  betafun <- discrete2curve(fdasrvf::beta[, , 1, 1])
  ctr <- get_curve_centroid(betafun)
  expect_length(ctr, nrow(fdasrvf::beta[, , 1, 1]))
  expect_true(all(is.finite(ctr)))
})
