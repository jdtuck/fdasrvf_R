q <- f_to_srvf(simu_data$f, simu_data$time)
time <- simu_data$time
pens <- c("none", "roughness", "l2gam", "l2psi", "geodesic")

test_that("`optimum.reparam()` returns a warping function for every penalty", {
  for (method in c("DP", "DPo", "RBFGS")) {
    for (pen in pens) {
      gam <- optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0.01,
                             pen = pen, method = method)
      expect_length(gam, length(time))
      expect_equal(gam[1], 0)
      expect_equal(gam[length(gam)], 1)
      expect_true(all(diff(gam) >= -1e-12))
    }
  }
})

test_that("`optimum.reparam()` penalty has no effect when lambda is 0", {
  for (method in c("DP", "DPo", "RBFGS")) {
    gam0 <- optimum.reparam(q[, 1], time, q[, 2], time, method = method)
    for (pen in pens) {
      gam <- optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0,
                             pen = pen, method = method)
      expect_equal(gam, gam0)
    }
  }
})

test_that("`optimum.reparam()` penalty 'none' ignores lambda", {
  for (method in c("DP", "DPo", "RBFGS")) {
    gam0 <- optimum.reparam(q[, 1], time, q[, 2], time, method = method)
    gam <- optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0.5,
                           pen = "none", method = method)
    expect_equal(gam, gam0)
  }
})

test_that("`optimum.reparam()` passes the penalty to the DP method", {
  gam <- lapply(c("none", "roughness", "l2gam"), function(pen) {
    optimum.reparam(q[, 1], time, q[, 2], time, lambda = 0.01, pen = pen,
                    method = "DP")
  })
  expect_false(isTRUE(all.equal(gam[[1]], gam[[2]])))
  expect_false(isTRUE(all.equal(gam[[2]], gam[[3]])))
})

test_that("`optimum.reparam()` rejects an unknown penalty", {
  expect_error(
    optimum.reparam(q[, 1], time, q[, 2], time, pen = "bogus"),
    "invalid penalty selection"
  )
})
