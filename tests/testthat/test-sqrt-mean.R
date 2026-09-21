# Tests for the Karcher-mean / median of warping functions computed on the
# Hilbert sphere of square-root densities: SqrtMean(), SqrtMedian() and
# SqrtMeanInverse() from R/SqrtMean.R, R/SqrtMedian.R and R/geometry.R.
# These are iterative and depend only on the (deterministic) input warping
# functions, so the assertions cover the return structure plus the invariant
# that the recovered mean/median warping function is itself a valid warping
# function (starts at 0, ends at 1, monotone non-decreasing).

is_valid_warping <- function(g, tol = 1e-6) {
  abs(g[1]) < 1e-3 &&
    abs(g[length(g)] - 1) < 1e-3 &&
    all(diff(g) >= -tol)
}

test_that("`SqrtMean()` works", {
  w <- simu_warp$warping_functions
  out <- SqrtMean(w)
  expect_equal(names(out), c("mu", "gam_mu", "psi", "vec"))
  expect_length(out$mu, nrow(w))
  expect_equal(dim(out$psi), dim(w))
  expect_equal(dim(out$vec), dim(w))
  expect_length(out$gam_mu, nrow(w))
  expect_true(is_valid_warping(out$gam_mu))
})

test_that("`SqrtMedian()` works", {
  w <- simu_warp$warping_functions
  out <- SqrtMedian(w)
  expect_equal(names(out), c("median", "gam_median", "psi", "vec"))
  expect_length(out$median, nrow(w))
  expect_length(out$gam_median, nrow(w))
  expect_true(is_valid_warping(out$gam_median))
})

test_that("`SqrtMeanInverse()` returns a valid warping function", {
  w <- simu_warp$warping_functions
  gamI <- SqrtMeanInverse(w)
  expect_length(gamI, nrow(w))
  expect_true(is_valid_warping(gamI))
})
