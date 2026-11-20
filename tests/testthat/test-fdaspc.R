skip_on_cran()

warp <- time_warping(simu_data$f, simu_data$time, lambda = 0)
t <- simu_data$time

test_that("in-sample T2 matches the scaled fPCA coefficients", {
  s <- fdaspc(warp, realign = FALSE)
  coef_amp <- matrix(s$vpca$coef, nrow = ncol(warp$f0))
  coef_ph <- matrix(s$hpca$coef, nrow = ncol(warp$f0))
  expect_equal(s$amp$T2,
               rowSums(sweep(coef_amp ^ 2, 2, s$vpca$latent, "/")),
               tolerance = 1e-10)
  expect_equal(s$ph$T2,
               rowSums(sweep(coef_ph ^ 2, 2, s$hpca$latent, "/")),
               tolerance = 1e-10)
})

test_that("contributions sum to T2 and SPE", {
  s <- fdaspc(warp, no = c(3, 2))
  expect_equal(c(s$amp$k, s$ph$k), c(3, 2))
  s <- spc_monitor(s, warp$f0[, 1:4])
  for (m in list(s$mon$amp, s$mon$ph)) {
    expect_equal(rowSums(m$contT2), m$T2, tolerance = 1e-8)
    expect_equal(rowSums(m$contSPE), m$SPE, tolerance = 1e-8)
  }
})

test_that("training data is mostly in control", {
  s <- spc_monitor(fdaspc(warp, alpha = 0.05), warp$f0)
  expect_lt(mean(s$mon$alarm), 0.5)
})

test_that("amplitude and phase faults raise alarms", {
  s <- fdaspc(warp, alpha = 0.05)
  f <- warp$f0[, 1:3]
  g <- pmin(pmax(t + 0.08 * sin(2 * pi * t), 0), 1)
  fp <- apply(f, 2, function(fi) warp_f_gamma(fi, t, g))
  s <- spc_monitor(s, cbind(1.5 * f, fp))
  expect_true(all(s$mon$amp$alarm[1:3]))
  expect_true(all(s$mon$ph$alarm[4:6]))
  expect_equal(s$mon$alarm, s$mon$amp$alarm | s$mon$ph$alarm)
})
