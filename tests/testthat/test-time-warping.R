test_that("The function `time_warping()` works", {
  out <- time_warping(
    f = simu_data$f,
    time = simu_data$time,
    max_iter = 1
  )
  expect_equal(length(out), 15)
  expect_equal(names(out), c("time", "f0", "q0", "fn", "qn", "fmean", "mqn",
                             "warping_functions", "original_variance",
                             "amplitude_variance", "phase_variance", "qun",
                             "inverse_average_warping_function", "rsamps",
                             "call"))
  expect_equal(length(out$time), 101)
  expect_equal(dim(out$f0), c(101, 21))
  expect_equal(dim(out$q0), c(101, 21))
  expect_equal(dim(out$fn), c(101, 21))
  expect_equal(dim(out$qn), c(101, 21))
  expect_equal(length(out$fmean), 101)
  expect_equal(length(out$mqn), 101)
  expect_equal(dim(out$warping_functions), c(101, 21))
  expect_equal(length(out$original_variance), 1)
  expect_equal(length(out$amplitude_variance), 1)
  expect_equal(length(out$phase_variance), 1)
  expect_equal(length(out$qun), 2)
  expect_equal(length(out$inverse_average_warping_function), 101)
  expect_equal(length(out$call), 9)
  expect_equal(names(out$call), c("lambda", "penalty_method", "centroid_type",
                                  "center_warpings", "smooth_data", "sparam",
                                  "parallel", "optim_method", "max_iter"))
  expect_snapshot(out)
})
