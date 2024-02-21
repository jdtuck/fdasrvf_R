test_that("`calc_shape_dist()` works", {
  beta1 <- fdasrvf::beta[, , 1, 1]
  R <- matrix(
    c(sqrt(2) / 2, -sqrt(2) / 2, sqrt(2) / 2, sqrt(2) / 2),
    nrow = 2, ncol = 2
  )
  beta1r <- R %*% beta1
  out <- calc_shape_dist(
    beta1 = beta1,
    beta2 = beta1r,
    mode = "O",
    alignment = FALSE,
    rotation = FALSE,
    scale = FALSE,
    lambda = 0
  )
  expect_true(out$d > 0)
  out <- calc_shape_dist(
    beta1 = beta1,
    beta2 = beta1r,
    mode = "O",
    alignment = FALSE,
    rotation = TRUE,
    scale = FALSE,
    lambda = 0
  )
  expect_equal(out$d, 0, tolerance = 1e-4)
  out <- calc_shape_dist(
    beta1 = beta1,
    beta2 = beta1r,
    mode = "O",
    alignment = TRUE,
    rotation = TRUE,
    scale = FALSE,
    lambda = 0
  )
  expect_equal(out$d, 0, tolerance = 1e-4)
  out <- calc_shape_dist(
    beta1 = beta1,
    beta2 = beta1r,
    mode = "O",
    alignment = FALSE,
    rotation = TRUE,
    scale = TRUE,
    lambda = 0
  )
  expect_equal(out$d, 0, tolerance = 1e-4)
})
