# Warping helpers on small grids. `smooth_warp()` is a diffeomorphism of the
# unit square onto itself (positive Jacobian, identity on the boundary).

smooth_warp <- function(m, n) {
  gam <- makediffeoid(m, n)
  X <- gam[, , 1]
  Y <- gam[, , 2]
  gam[, , 1] <- X + 0.1 * sin(pi * X) * sin(pi * Y)
  gam[, , 2] <- Y + 0.05 * sin(2 * pi * X) * sin(pi * Y)
  gam
}

grid_sizes <- list(c(12, 12), c(8, 13), c(13, 8))

test_that("`makediffeoid()` puts column coordinates first", {
  gam <- makediffeoid(3, 5)
  expect_equal(dim(gam), c(3, 5, 2))
  expect_equal(gam[1, , 1], seq(0, 1, length.out = 5))
  expect_equal(gam[, 1, 2], seq(0, 1, length.out = 3))
})

test_that("the identity warp reproduces the image", {
  set.seed(1)
  for (sz in grid_sizes) {
    id <- makediffeoid(sz[1], sz[2])
    img <- matrix(runif(prod(sz)), sz[1], sz[2])
    expect_equal(apply_gam_to_imag(img, id), img)
    img3 <- array(runif(prod(sz) * 3), c(sz, 3))
    expect_equal(apply_gam_to_imag(img3, id), img3)
  }
})

test_that("composing with the identity leaves a warp unchanged", {
  for (sz in grid_sizes) {
    id <- makediffeoid(sz[1], sz[2])
    gam <- smooth_warp(sz[1], sz[2])
    expect_equal(apply_gam_to_gam(id, id), id)
    expect_equal(apply_gam_to_gam(id, gam), gam)
    expect_equal(apply_gam_to_gam(gam, id), gam)
    expect_equal(apply_gam_gamid(id, gam), gam)
  }
})

test_that("warping evaluates the image at the paired warp coordinates", {
  # The spline reproduces a plane exactly, so img o gam is known in closed
  # form.
  for (sz in grid_sizes) {
    id <- makediffeoid(sz[1], sz[2])
    plane <- function(g) 1 + 2 * g[, , 1] - 3 * g[, , 2]
    gam <- smooth_warp(sz[1], sz[2])
    expect_equal(apply_gam_to_imag(plane(id), gam), plane(gam))
  }
})

test_that("warping clamps points to the unit square", {
  id <- makediffeoid(6, 7)
  img <- 1 + 2 * id[, , 1] - 3 * id[, , 2]
  gam <- id
  gam[, , 1] <- gam[, , 1] + 0.2
  out <- apply_gam_to_imag(img, gam)
  expect_false(anyNA(out))
  expect_equal(out[, 7], img[, 7])
})

test_that("warping can resample onto a grid of a different size", {
  id_small <- makediffeoid(5, 6)
  id_big <- makediffeoid(9, 11)
  img <- 1 + 2 * id_big[, , 1] - 3 * id_big[, , 2]
  expect_equal(apply_gam_to_imag(img, id_small),
               1 + 2 * id_small[, , 1] - 3 * id_small[, , 2])
})

test_that("`gradient2()` handles non-square images", {
  id <- makediffeoid(6, 9)
  img <- 2 * id[, , 1] - 3 * id[, , 2]
  out <- gradient2(img, 1 / 8, 1 / 5)
  expect_equal(out$dxdu, matrix(2, 6, 9))
  expect_equal(out$dydv, matrix(-3, 6, 9))
})

test_that("`pair_align_image()` runs on square and non-square images", {
  for (rows in list(1:64, 9:40)) {
    I1 <- im$I1[rows, 5:60]
    I2 <- im$I2[rows, 5:60]
    expect_output(out <- pair_align_image(I1, I2, itermax = 2), "Iteration 0")
    expect_equal(dim(out$I2_new), dim(I1))
    expect_equal(dim(out$gam), c(dim(I1), 2))
    expect_false(anyNA(out$I2_new))
    expect_true(check_crossing(out$gam))
  }
})
