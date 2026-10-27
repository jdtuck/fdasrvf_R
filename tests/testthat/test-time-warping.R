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
  expect_equal(length(out$qun), 1)
  expect_equal(length(out$inverse_average_warping_function), 101)
  expect_equal(length(out$call), 9)
  expect_equal(names(out$call), c("lambda", "penalty_method", "centroid_type",
                                  "center_warpings", "smooth_data", "sparam",
                                  "parallel", "optim_method", "max_iter"))
  expect_snapshot(out)
})

warp_by <- function(f, gam, time) {
  M <- length(time)
  sapply(seq_len(ncol(f)), function(n)
    stats::approx(time, f[, n], xout = (time[M] - time[1]) * gam[, n] + time[1])$y)
}

test_that("`time_warping()` returns fn = f o warping_functions, as MATLAB does", {
  f <- simu_data$f[, 1:8]
  time <- simu_data$time
  for (args in list(list(), list(centroid_type = "median"),
                    list(center_warpings = FALSE))) {
    out <- do.call(time_warping, c(list(f = f, time = time), args))
    expect_equal(out$fn, warp_by(f, out$warping_functions, time),
                 tolerance = 1e-12)
    expect_equal(out$qn, f_to_srvf(out$fn, time), tolerance = 1e-12)
    if (length(args) == 0)
      expect_equal(out$mqn, rowMeans(out$qn), tolerance = 1e-12)
    # qun holds one relative template change per iteration
    expect_true(all(out$qun >= 0))
    expect_lte(length(out$qun), out$call$max_iter)
  }
})

test_that("`time_warping()` matches the MATLAB reference statistics", {
  # MATLAB fdawarp(f, t).time_warping(0, method = 'DP') on simu_data gave
  # amp_var 0.0200 and qun = [0.0551 0.0099] with the piecewise-constant DP
  # cost; with the cost integrating the SRVFs linearly and nbhd_dim = 10,
  # R and Python fdasrsf (srsf_align) give amp_var 0.0200446 and
  # qun = [0.0549 0.0091]
  out <- time_warping(simu_data$f, simu_data$time,
                      penalty_method = "l2psi")
  expect_equal(out$amplitude_variance, 0.0200446, tolerance = 1e-5)
  expect_equal(out$qun, c(0.0549, 0.0091), tolerance = 1e-2)
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
