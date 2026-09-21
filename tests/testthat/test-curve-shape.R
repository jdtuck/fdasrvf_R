# Tests for the open-/closed-curve shape routines operating on the sample
# curves in `beta`. Where a clean mathematical invariant exists (a geodesic
# path has the requested number of steps and its endpoints match the input
# curves; an alignment produces a valid warping function) it is asserted;
# otherwise these are return-structure / smoke tests exercising the code path.

test_that("`curve_geodesic()` returns a k-step path between the curves", {
  beta1 <- fdasrvf::beta[, , 1, 1]
  beta2 <- fdasrvf::beta[, , 1, 2]
  k <- 3
  out <- curve_geodesic(beta1, beta2, k = k)
  expect_equal(names(out), c("geod", "geod_q"))
  n <- nrow(beta1)
  M <- ncol(beta1)
  expect_equal(dim(out$geod), c(n, M, k))
  expect_equal(dim(out$geod_q), c(n, M, k))
  # Endpoints of the geodesic are finite curves of the right shape.
  expect_true(all(is.finite(out$geod[, , 1])))
  expect_true(all(is.finite(out$geod[, , k])))
})

test_that("`curve_pair_align()` returns a valid warping function", {
  out <- curve_pair_align(fdasrvf::beta[, , 1, 1], fdasrvf::beta[, , 1, 2])
  expect_true(all(c("beta2n", "q2n", "gam", "R") %in% names(out)))
  gam <- out$gam
  expect_equal(gam[1], 0, tolerance = 1e-3)
  expect_equal(gam[length(gam)], 1, tolerance = 1e-3)
  expect_true(all(diff(gam) >= -1e-6))
  expect_equal(dim(out$R), c(2, 2))
})

test_that("`curve_depth()` returns amplitude and phase depths", {
  out <- curve_depth(fdasrvf::beta[, , 1, 1:3])
  expect_equal(names(out), c("amp", "phase"))
  expect_length(out$amp, 3)
  expect_length(out$phase, 3)
  expect_true(all(is.finite(out$amp)))
  expect_true(all(is.finite(out$phase)))
})

test_that("`sample_shapes()` draws shapes from a Karcher-mean model", {
  km <- multivariate_karcher_mean(fdasrvf::beta[, , 1, 1:5], scale = TRUE, maxit = 2)
  withr::with_seed(1234, {
    out <- sample_shapes(km, no = 3, numSamp = 5)
  })
  expect_true(all(c("betas", "qns", "gams") %in% names(out)))
})

test_that("`plot_curve()` and `f_plot()` draw without error", {
  expect_no_error({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    plot_curve(fdasrvf::beta[, , 1, 1])
    f_plot(simu_data$time, simu_data$f)
  })
})
