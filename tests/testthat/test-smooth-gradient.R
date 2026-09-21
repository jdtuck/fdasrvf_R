# Tests for the low-level numeric utilities:
#   * gradient()      (R/gradient.R)   - finite-difference derivative,
#   * smooth.data()   (R/smooth.data.R) - moving-average smoother,
#   * resamplecurve() (R/resamplecurve.R) - arc-length resampling,
#   * rgam()          (R/rgam.R)        - random warping-function generator,
#   * interparc()     (R/interparc.R)   - arc-length interpolation.
# Correctness is anchored on analytic derivatives (the gradient of a linear
# ramp is its slope; the gradient of a constant is zero), shape preservation,
# and the validity of generated warping functions.

test_that("`gradient()` matches known analytic derivatives", {
  M <- 101
  time <- seq(0, 1, length.out = M)
  binsize <- mean(diff(time))

  # d/dt of the identity ramp is 1 everywhere.
  expect_equal(gradient(time, binsize), rep(1, M), tolerance = 1e-8)
  # d/dt of a constant is 0 everywhere.
  expect_equal(gradient(rep(3.5, M), binsize), rep(0, M), tolerance = 1e-8)
  # Length is preserved.
  expect_length(gradient(sin(2 * pi * time), binsize), M)
})

test_that("`gradient()` preserves matrix shape for multidimensional input", {
  x <- fdasrvf::beta[, , 1, 1]
  binsize <- 1 / (ncol(x) - 1)
  g <- gradient(x, binsize, multidimensional = TRUE)
  expect_equal(dim(g), dim(x))
})

test_that("`smooth.data()` preserves shape and stays finite", {
  f <- simu_data$f
  out <- smooth.data(f, sparam = 5)
  expect_equal(dim(out), dim(f))
  expect_true(all(is.finite(out)))
})

test_that("`resamplecurve()` resamples to the requested number of points", {
  x <- fdasrvf::beta[, , 1, 1]
  out <- resamplecurve(x, N = 50)
  expect_equal(dim(out), c(nrow(x), 50))
  expect_true(all(is.finite(out)))
})

test_that("`rgam()` generates valid random warping functions", {
  withr::with_seed(1234, {
    gam <- rgam(N = 101, sigma = 0.1, num = 3)
  })
  expect_equal(dim(gam), c(3, 101))
  for (k in 1:3) {
    g <- gam[k, ]
    expect_equal(g[1], 0, tolerance = 1e-6)
    expect_equal(g[length(g)], 1, tolerance = 1e-6)
    expect_true(all(diff(g) >= -1e-6))
  }
})

test_that("`interparc()` interpolates along arc length", {
  x <- cos(seq(0, 2 * pi, length.out = 20) * 1.5 - pi / 2)
  y <- sin(seq(0, 2 * pi, length.out = 20))
  out <- interparc(10, x, y)
  expect_equal(nrow(out), 10)
  expect_true(all(is.finite(as.matrix(out))))
})
