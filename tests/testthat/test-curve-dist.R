test_that("`curve_dist()` works", {
  N <- 4
  betas <- fdasrvf::beta[, , 1, 1:N]
  out <- curve_dist(
    beta = betas,
    mode = "O",
    alignment = FALSE,
    rotation = FALSE,
    scale = FALSE,
    include.length = FALSE,
    lambda = 0.0
  )
  expect_equal(attr(out$Da, "Size"), N)
  out <- curve_dist(
    beta = betas,
    mode = "O",
    alignment = TRUE,
    rotation = FALSE,
    scale = FALSE,
    include.length = FALSE,
    lambda = 0.0
  )
  expect_equal(attr(out$Da, "Size"), N)
  out <- curve_dist(
    beta = betas,
    mode = "O",
    alignment = FALSE,
    rotation = TRUE,
    scale = FALSE,
    include.length = FALSE,
    lambda = 0.0
  )
  expect_equal(attr(out$Da, "Size"), N)
  out <- curve_dist(
    beta = betas,
    mode = "O",
    alignment = FALSE,
    rotation = FALSE,
    scale = TRUE,
    include.length = FALSE,
    lambda = 0.0
  )
  expect_equal(attr(out$Da, "Size"), N)
  out <- curve_dist(
    beta = betas,
    mode = "O",
    alignment = FALSE,
    rotation = FALSE,
    scale = TRUE,
    include.length = TRUE,
    lambda = 0.0
  )
  expect_equal(attr(out$Da, "Size"), N)
})
