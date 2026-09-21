# Tests for the discrete warping-application helpers warp_f_gamma() and
# warp_q_gamma(). The key invariant is that warping by the identity warping
# gamma(t) = t leaves the function (resp. SRVF) unchanged, and that the output
# has the same length as the input for both interpolation methods.

test_that("`warp_f_gamma()` under the identity warping is a no-op", {
  f <- simu_data$f[, 1]
  time <- simu_data$time
  gamid <- seq(0, 1, length.out = length(time))

  f_lin <- warp_f_gamma(f, time, gamid, spl.int = FALSE)
  f_spl <- warp_f_gamma(f, time, gamid, spl.int = TRUE)

  expect_length(f_lin, length(f))
  expect_length(f_spl, length(f))
  expect_equal(f_lin, f, tolerance = 1e-6)
  expect_equal(f_spl, f, tolerance = 1e-4)
})

test_that("`warp_q_gamma()` under the identity warping is a no-op", {
  f <- simu_data$f[, 1]
  time <- simu_data$time
  q <- f_to_srvf(f, time)
  gamid <- seq(0, 1, length.out = length(time))

  q_lin <- warp_q_gamma(q, time, gamid, spl.int = FALSE)
  q_spl <- warp_q_gamma(q, time, gamid, spl.int = TRUE)

  expect_length(q_lin, length(q))
  expect_length(q_spl, length(q))
  # sqrt(gamma') = 1 for the identity, so the SRVF is unchanged.
  expect_equal(q_lin, q, tolerance = 1e-6)
  expect_equal(q_spl, q, tolerance = 1e-4)
})
