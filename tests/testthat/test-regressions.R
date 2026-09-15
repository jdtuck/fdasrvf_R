test_that("interp1_flat interpolates between two flat segments", {
  x <- c(0, 0.2, 0.4, 0.4, 0.6, 0.8, 0.8, 1)
  y <- c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 1)
  expect_equal(interp1_flat(x, y, c(0.5, 0.7)), c(0.4, 0.6))
})

test_that("simpson integrates every column of a matrix", {
  x <- seq(0, 1, length.out = 10)
  y <- cbind(x^2, 3 * x^2)
  expect_equal(simpson(x, y), c(1 / 3, 1))
  expect_equal(simpson(x, y), c(simpson(x, y[, 1]), simpson(x, y[, 2])))
})

test_that("gradient2 handles non-square matrices", {
  a <- outer(1:5, 1:3)
  out <- gradient2(a)
  expect_equal(dim(out$dydv), c(5, 3))
  expect_equal(out$dydv[, 2], rep(2, 5))
  expect_equal(out$dxdu[3, ], rep(3, 3))
})

test_that("optimum.reparam only uses DPo when the grids match", {
  q <- f_to_srvf(simu_data$f[, 1:2], simu_data$time)
  t2 <- simu_data$time
  t2[50] <- t2[50] + 1e-4
  gam <- optimum.reparam(q[, 1], simu_data$time, q[, 2], t2, method = "DPo")
  expect_equal(gam, optimum.reparam(q[, 1], simu_data$time, q[, 2], t2,
                                    method = "DP"))
})

test_that("get_distance_matrix honours scale = TRUE", {
  q1 <- curve2srvf(fdasrvf::beta[, , 1, 1])
  q2 <- curve2srvf(fdasrvf::beta[, , 1, 2])
  dm <- get_distance_matrix(list(q1, q2), scale = TRUE)
  expect_equal(
    as.numeric(dm$Da),
    get_shape_distance(q1, q2, scale = TRUE)$amplitude_distance,
    tolerance = 1e-6
  )
})

test_that("multivariate_karcher_mean handles convergence on the first iteration", {
  out <- multivariate_karcher_mean(fdasrvf::beta[, , 1, c(1, 1, 1)], maxit = 5)
  expect_equal(dim(out$v), dim(out$q))
  expect_equal(out$qn, out$q)
})

test_that("predict.curve_pca reproduces the coefficients of the fitted curves", {
  km <- multivariate_karcher_mean(fdasrvf::beta[, , 1, 1:4], maxit = 20)
  pc <- multivariate_pca(km, no = 2, showplot = FALSE)
  a <- predict(pc)
  expect_equal(dim(a), c(2, 4))
  expect_equal(a, pc$coef, tolerance = 1e-6)
  expect_equal(predict(pc, fdasrvf::beta[, , 1, 2]), pc$coef[, 2, drop = FALSE],
               tolerance = 1e-6)
})

test_that("v_to_curve maps each column of a matrix of shooting vectors", {
  mu <- curve_to_q(fdasrvf::beta[, , 1, 1])$q
  V <- matrix(0, length(mu), 3)
  out <- v_to_curve(V, mu)
  expect_equal(dim(out), dim(V))
  expect_equal(out[, 2], c(q_to_curve(mu)))
})

test_that("getPersistentPeaks clusters peaks by how often they persist", {
  counts <- rbind(rep(1, 5), rep(1, 5), c(1, NaN, NaN, NaN, NaN), rep(1, 5))
  expect_equal(getPersistentPeaks(counts), c(1L, 2L, 4L))
})

test_that("gauss_model returns M x n warps and supports sort_samples", {
  M <- length(simu_warp$time)
  set.seed(1)
  out <- gauss_model(simu_warp, n = 3, sort_samples = TRUE)
  expect_equal(dim(out$ft), c(M, 3))
  expect_equal(dim(out$gams), c(M, 3))
  expect_false(anyNA(out$ft))
})

test_that("reparam_curve with DPo returns the optimal rotation", {
  out <- reparam_curve(fdasrvf::beta[, , 1, 1], fdasrvf::beta[, , 1, 2], method = "DPo")
  expect_equal(dim(out$R), c(2, 2))
})

test_that("kmeans_align copes with a single-curve cluster", {
  f <- cbind(simu_data$f[, 1:5], 10 * simu_data$f[, 6])
  out <- kmeans_align(f, simu_data$time, K = 2, seeds = c(1, 6), max_iter = 2)
  expect_equal(out$labels[6], 2)
  expect_equal(sum(out$labels == 2), 1)
})

test_that("kmeans_align keeps multivariate medoid templates in the right slots", {
  out <- kmeans_align(fdasrvf::beta[, , 1, 1:4], seq(0, 1, length.out = 100), K = 1,
                      seeds = 1, centroid_type = "medoid", scale = FALSE,
                      max_iter = 1)
  expect_equal(out$templates.q[, , 1],
               curve_to_q(out$templates[, , 1], scale = FALSE)$q)
})

test_that("jacob_imag uses all three cross-product terms for 3-d maps", {
  g <- makediffeoid(8, 8)
  F1 <- array(0, c(8, 8, 3))
  F1[, , 1] <- g[, , 1]
  F1[, , 3] <- g[, , 2]
  expect_true(all(jacob_imag(F1) > 0))
})

test_that("the 's' image basis varies with its own frequency and axis", {
  b <- formbasisTid(2, 8, 8, "s")$b
  # the second component of the first element depends on the row only
  expect_equal(b[, 1, 2, 1], b[, 5, 2, 1])
  expect_lt(abs(stats::cor(c(b[, , 2, 1]), c(b[, , 2, 3]))), 0.99)
})

test_that("time_warping keeps iterating until the template settles", {
  out <- time_warping(simu_data$f, simu_data$time)
  n <- length(out$qun)
  expect_true(out$qun[n - 1] < 1e-2 || n == 20 + 2)
})

test_that("calc_j fills in the (2, 2) entry", {
  b <- list(matrix(1, 2, 5), matrix(2, 2, 5))
  j <- calc_j(b)
  expect_equal(j[2, 2], 8)
  expect_equal(j[1, 2], 4)
})
