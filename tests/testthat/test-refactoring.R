test_that("discrete to continuous functions work", {
  N <- dim(fdasrvf::beta)[4]
  beta_mats <- fdasrvf::beta[, , 1, ]
  gamma_vec <- toy_warp$gam[, 1]

  betafuns <- lapply(1:N, \(n) discrete2curve(beta_mats[, , n]))
  expect_equal(dim(betafuns[[1]](c(0.2, 0.4, 0.6), deriv = 1)), c(2, 3))
  expect_equal(dim(betafuns[[1]](0.4)), c(2, 1))

  gamfun <- discrete2warping(gamma_vec)
  expect_equal(length(gamfun(c(0.2, 0.4, 0.6))), 3)
  expect_equal(length(gamfun(0.4)), 1)
})

test_that("`curve2srvf()` works", {
  N <- dim(fdasrvf::beta)[4]
  beta_mats <- fdasrvf::beta[, , 1, ]
  gamma_vec <- toy_warp$gam[, 1]

  betafuns <- lapply(1:N, \(n) discrete2curve(beta_mats[, , n]))
  qfuns <- lapply(betafuns, curve2srvf)
  expect_equal(dim(qfuns[[1]](c(0.2, 0.4, 0.6))), c(2, 3))
  expect_equal(dim(qfuns[[1]](0.4)), c(2, 1))
})

test_that("`curve2srvf()` and `srvf2curve()` are inverses", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  idx1 <- 1
  idx2 <- 6

  q <- curve2srvf(fdasrvf::beta[, , 1, 1])
  beta_recon <- srvf2curve(q, beta0 = fdasrvf::beta[, 1, 1, 1])
  expect_true(get_l2_distance(betafuns[[1]], beta_recon) < 6e-6)

  q <- curve2srvf(betafuns[[1]])
  beta_recon <- srvf2curve(q, beta0 = fdasrvf::beta[, 1, 1, 1])
  expect_true(get_l2_distance(betafuns[[1]], beta_recon) < 6e-6)
})

test_that("`get_l2_distance()` is symmetric", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6

  d1 <- get_l2_distance(qfuns[[idx1]], qfuns[[idx2]])
  d2 <- get_l2_distance(qfuns[[idx2]], qfuns[[idx1]])
  expect_true(abs(d1 - d2) < .Machine$double.eps)
})

test_that("`warp_srvf()` works", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6
  gamfun <- discrete2warping(toy_warp$gam[, 1])

  q1 <- warp_srvf(
    qfun = qfuns[[idx1]],
    gamfun = gamfun,
    betafun = betafuns[[idx1]]
  )
  q2 <- warp_srvf(
    qfun = qfuns[[idx1]],
    gamfun = gamfun
  )
  expect_true(get_l2_distance(q1, q2) < 5e-15)
})

test_that("`get_l2_distance()` is Gamma-invartiant", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6
  gamfun <- discrete2warping(toy_warp$gam[, 1])

  q1 <- warp_srvf(
    qfun = qfuns[[idx1]],
    gamfun = gamfun,
    betafun = betafuns[[idx1]]
  )
  q2 <- warp_srvf(
    qfun = qfuns[[idx2]],
    gamfun = gamfun,
    betafun = betafuns[[idx2]]
  )

  d1 <- get_l2_distance(qfuns[[idx1]], qfuns[[idx2]])
  d2 <- get_l2_distance(q1, q2)
  expect_true(abs(d1 - d2) < 3e-7)
})

test_that("`get_l2_inner_product()` is symmetric", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6

  ip1 <- get_l2_inner_product(qfuns[[idx1]], qfuns[[idx2]])
  ip2 <- get_l2_inner_product(qfuns[[idx2]], qfuns[[idx1]])
  expect_true(abs(ip1 - ip2) < 1e-4)
})

test_that("`get_l2_inner_product()` is Gamma-invariant", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6
  gamfun <- discrete2warping(toy_warp$gam[, 1])

  q1 <- warp_srvf(
    qfun = qfuns[[idx1]],
    gamfun = gamfun,
    betafun = betafuns[[idx1]]
  )
  q2 <- warp_srvf(
    qfun = qfuns[[idx2]],
    gamfun = gamfun,
    betafun = betafuns[[idx2]]
  )

  ip1 <- get_l2_inner_product(qfuns[[idx1]], qfuns[[idx2]])
  ip2 <- get_l2_inner_product(q1, q2)
  expect_true(abs(ip1 - ip2) < 7e-4)
})

test_that("`get_warping_distance()` is symmetric", {
  gam1 <- discrete2warping(toy_warp$gam[, 1])
  gam2 <- discrete2warping(toy_warp$gam[, 2])

  d1 <- get_warping_distance(gam1, gam2)
  d2 <- get_warping_distance(gam2, gam1)
  expect_true(abs(d1 - d2) < 9e-7)

  d1 <- get_warping_distance(get_identity_warping(), gam1)
  d2 <- get_warping_distance(gam1, get_identity_warping())
  expect_true(abs(d1 - d2) < 2e-6)
})

test_that("`get_hilbert_sphere_distance()` is symmetric", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6

  q1p <- to_hilbert_sphere(qfuns[[idx1]])
  q2p <- to_hilbert_sphere(qfuns[[idx2]])

  d1 <- get_hilbert_sphere_distance(q1p, q2p)
  d2 <- get_hilbert_sphere_distance(q2p, q1p)
  expect_true(abs(d1 - d2) < .Machine$double.eps)
})

test_that("`get_hilbert_sphere_distance()` is Gamma-invariant", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6
  gamfun <- discrete2warping(toy_warp$gam[, 1])

  q1p <- to_hilbert_sphere(qfuns[[idx1]])
  q2p <- to_hilbert_sphere(qfuns[[idx2]])

  d1 <- get_hilbert_sphere_distance(q1p, q2p)

  q1 <- warp_srvf(
    qfun = qfuns[[idx1]],
    gamfun = gamfun,
    betafun = betafuns[[idx1]]
  )
  q2 <- warp_srvf(
    qfun = qfuns[[idx2]],
    gamfun = gamfun,
    betafun = betafuns[[idx2]]
  )

  q1p <- to_hilbert_sphere(q1)
  q2p <- to_hilbert_sphere(q2)
  d2 <- get_hilbert_sphere_distance(q1p, q2p)

  expect_true(abs(d1 - d2) < 4e-7)
})

test_that("`get_shape_distance()` is symmetric", {
  N <- dim(fdasrvf::beta)[4]
  betafuns <- lapply(1:N, \(n) discrete2curve(fdasrvf::beta[, , 1, n]))
  qfuns <- lapply(betafuns, curve2srvf)
  idx1 <- 1
  idx2 <- 6

  d1 <- get_shape_distance(qfuns[[idx1]], qfuns[[idx2]])
  d1 <- c(d1$amplitude_distance, d1$phase_distance)
  d2 <- get_shape_distance(qfuns[[idx2]], qfuns[[idx1]])
  d2 <- c(d2$amplitude_distance, d2$phase_distance)
  expect_true(all(abs(d1 - d2) < .Machine$double.eps))

  d1 <- get_shape_distance(qfuns[[idx1]], qfuns[[idx2]], alignment = TRUE)
  d1 <- c(d1$amplitude_distance, d1$phase_distance)
  d2 <- get_shape_distance(qfuns[[idx2]], qfuns[[idx1]], alignment = TRUE)
  d2 <- c(d2$amplitude_distance, d2$phase_distance)
  expect_true(all(abs(d1 - d2) < 3e-6))

  d1 <- get_shape_distance(qfuns[[idx1]], qfuns[[idx2]], rotation = TRUE)
  d1 <- c(d1$amplitude_distance, d1$phase_distance)
  d2 <- get_shape_distance(qfuns[[idx2]], qfuns[[idx1]], rotation = TRUE)
  d2 <- c(d2$amplitude_distance, d2$phase_distance)
  expect_true(all(abs(d1 - d2) < 6e-14))

  d1 <- get_shape_distance(qfuns[[idx1]], qfuns[[idx2]], scale = TRUE)
  d1 <- c(d1$amplitude_distance, d1$phase_distance)
  d2 <- get_shape_distance(qfuns[[idx2]], qfuns[[idx1]], scale = TRUE)
  d2 <- c(d2$amplitude_distance, d2$phase_distance)
  expect_true(all(abs(d1 - d2) < 1e-6))
})

test_that("`get_distance_matrix()` works", {
  N <- 4L
  srvfs <- lapply(1:N, \(n) curve2srvf(fdasrvf::beta[, , 1, n]))
  out <- get_distance_matrix(srvfs)
  expect_equal(class(out), "list")
  expect_equal(length(out), 2L)
  expect_equal(names(out), c("Da", "Dp"))
  expect_equal(class(out$Da), "dist")
  expect_equal(class(out$Dp), "dist")
  expect_equal(attr(out$Da, "Size"), N)
})
