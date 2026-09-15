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
  expect_equal(length(out$qun), 3)
  expect_equal(length(out$inverse_average_warping_function), 101)
  expect_equal(length(out$call), 9)
  expect_equal(names(out$call), c("lambda", "penalty_method", "centroid_type",
                                  "center_warpings", "smooth_data", "sparam",
                                  "parallel", "optim_method", "max_iter"))
  expect_snapshot(out)
})

pens <- c("roughness", "l2gam", "l2psi", "geodesic", "none", "norm")

test_that("`time_warping()` runs with every penalty method", {
  for (pen in pens) {
    out <- time_warping(f = simu_data$f, time = simu_data$time, lambda = 0.01,
                        penalty_method = pen, max_iter = 1)
    expect_s3_class(out, "fdawarp")
    expect_equal(dim(out$warping_functions), c(101, 21))
    expect_true(all(is.finite(out$warping_functions)))
  }
})

test_that("`time_warping()` treats penalty 'norm' as 'l2gam'", {
  out_norm <- time_warping(f = simu_data$f, time = simu_data$time,
                           lambda = 0.01, penalty_method = "norm",
                           max_iter = 1)
  out_l2gam <- time_warping(f = simu_data$f, time = simu_data$time,
                            lambda = 0.01, penalty_method = "l2gam",
                            max_iter = 1)
  expect_equal(out_norm$call$penalty_method, "l2gam")
  expect_equal(out_norm$warping_functions, out_l2gam$warping_functions)
})

test_that("`time_warping()` and `ppd()` only offer methods `optimum.reparam()` accepts", {
  reparam_methods <- eval(formals(optimum.reparam)$method)
  expect_true(all(eval(formals(time_warping)$optim_method) %in% reparam_methods))
  expect_true(all(eval(formals(ppd)$optim_method) %in% reparam_methods))
})

test_that("every `time_warping()` optim_method runs and is stored for reuse", {
  for (m in eval(formals(time_warping)$optim_method)) {
    out <- time_warping(f = simu_data$f, time = simu_data$time,
                        optim_method = m, max_iter = 1)
    expect_equal(out$call$optim_method, m)
    expect_true(all(is.finite(out$warping_functions)))
    # predict methods pass the stored method straight to optimum.reparam()
    gam <- optimum.reparam(out$mqn, out$time, out$qn[, 1], out$time,
                           method = out$call$optim_method)
    expect_true(all(is.finite(gam)))
  }
})

test_that("`time_warping()` and `ppd()` reject the removed 'DP2' method up front", {
  expect_error(
    time_warping(f = simu_data$f, time = simu_data$time,
                 optim_method = "DP2", max_iter = 1),
    class = "rlang_error"
  )
  expect_error(
    ppd(simu_data$f, simu_data$time, parallel = FALSE, optim_method = "DP2"),
    class = "rlang_error"
  )
})
