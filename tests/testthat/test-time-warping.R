test_that("The function `time_warping()` works", {
  out <- time_warping(
    f = simu_data$f,
    time = simu_data$time,
    max_iter = 1
  )
  expect_equal(length(out), 17)
  expect_equal(names(out), c("f0", "time", "fn", "qn", "q0", "fmean", "mqn",
                             "gam", "orig.var", "amp.var", "phase.var", "qun",
                             "lambda", "method", "omethod", "gamI", "rsamps"))
  expect_equal(dim(out$f0), c(101, 21))
  expect_equal(length(out$time), 101)
  expect_equal(dim(out$fn), c(101, 21))
  expect_equal(dim(out$qn), c(101, 21))
  expect_equal(dim(out$q0), c(101, 21))
  expect_equal(length(out$fmean), 101)
  expect_equal(length(out$mqn), 101)
  expect_equal(dim(out$gam), c(101, 21))
  expect_equal(length(out$orig.var), 1)
  expect_equal(length(out$amp.var), 1)
  expect_equal(length(out$phase.var), 1)
  expect_equal(length(out$qun), 2)
  expect_equal(length(out$lambda), 1)
  expect_equal(length(out$method), 1)
  expect_equal(length(out$omethod), 1)
  expect_equal(length(out$gamI), 101)
  expect_equal(length(out$rsamps), 1)
  expect_snapshot(out)
})
