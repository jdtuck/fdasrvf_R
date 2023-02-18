test_that("`kmeans_align()` works", {
  withr::with_seed(1234, {
    out <- kmeans_align(simu_data$f, time = simu_data$time, K = 1)
  })
  expect_equal(
    names(out),
    c("f0", "q0", "time", "fn", "qn", "gam", "labels",
      "templates", "templates.q", "lambda", "omethod", "qun")
  )
  expect_equal(dim(out$f0), c(101, 21))
  expect_equal(dim(out$q0), c(101, 21))
  expect_equal(length(out$time), 101)
  expect_equal(length(out$fn), 1)
  expect_equal(dim(out$fn[[1]]), c(101, 21))
  expect_equal(length(out$qn), 1)
  expect_equal(dim(out$qn[[1]]), c(101, 21))
  expect_equal(length(out$gam), 1)
  expect_equal(dim(out$gam[[1]]), c(101, 21))
  expect_equal(out$labels, rep(1, 21))
  expect_equal(dim(out$templates), c(1, 101, 1))
  expect_equal(dim(out$templates.q), c(1, 101, 1))
  expect_equal(out$lambda, 0)
  expect_equal(out$omethod, "DP")
  expect_snapshot(out)
})
