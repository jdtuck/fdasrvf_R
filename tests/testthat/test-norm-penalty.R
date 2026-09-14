# "norm" is kept as an alias for the "l2gam" penalty. Use lambda > 0 so the
# penalty actually enters the alignment.

test_that("`pair_align_functions()` treats penalty 'norm' as 'l2gam'", {
  f1 <- simu_data$f[, 1]
  f2 <- simu_data$f[, 2]
  out_norm <- pair_align_functions(f1, f2, simu_data$time, lambda = 0.01,
                                   pen = "norm")
  out_l2gam <- pair_align_functions(f1, f2, simu_data$time, lambda = 0.01,
                                    pen = "l2gam")
  expect_true(all(is.finite(out_norm$gam)))
  expect_equal(out_norm, out_l2gam)
})

test_that("`multiple_align_functions()` treats penalty 'norm' as 'l2gam'", {
  f <- simu_data$f
  mu <- rowMeans(f)
  out_norm <- multiple_align_functions(f, simu_data$time, mu, lambda = 0.01,
                                       pen = "norm", showplot = FALSE,
                                       verbose = FALSE)
  out_l2gam <- multiple_align_functions(f, simu_data$time, mu, lambda = 0.01,
                                        pen = "l2gam", showplot = FALSE,
                                        verbose = FALSE)
  expect_s3_class(out_norm, "fdawarp")
  expect_equal(out_norm$call$penalty_method, "l2gam")
  expect_equal(out_norm$warping_functions, out_l2gam$warping_functions)
})

test_that("`elastic.distance()` treats penalty 'norm' as 'l2gam'", {
  f1 <- simu_data$f[, 1]
  f2 <- simu_data$f[, 2]
  out_norm <- elastic.distance(f1, f2, simu_data$time, lambda = 0.01,
                               pen = "norm")
  out_l2gam <- elastic.distance(f1, f2, simu_data$time, lambda = 0.01,
                                pen = "l2gam")
  expect_true(is.finite(out_norm$Dy) && is.finite(out_norm$Dx))
  expect_equal(out_norm, out_l2gam)
})

test_that("`elastic.depth()` treats penalty 'norm' as 'l2gam'", {
  f <- simu_data$f[, 1:4]
  out_norm <- elastic.depth(f, simu_data$time, lambda = 0.01, pen = "norm")
  out_l2gam <- elastic.depth(f, simu_data$time, lambda = 0.01, pen = "l2gam")
  expect_true(all(is.finite(out_norm$amp)) && all(is.finite(out_norm$phase)))
  expect_equal(out_norm, out_l2gam)
})
