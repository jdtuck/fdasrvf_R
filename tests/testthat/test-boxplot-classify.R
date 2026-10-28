out_median <- time_warping(simu_data$f, simu_data$time,
                           centroid_type = "median", max_iter = 1L)
amp_box <- boxplot(out_median, variability_type = "amplitude", what = "stats")
ph_box <- boxplot(out_median, variability_type = "phase", what = "stats")

test_that("boxplots store the region edges they use", {
  for (b in list(amp_box, ph_box)) {
    N <- length(b$dist)
    d <- sort(b$dist)
    expect_equal(b$alpha, 0.05)
    expect_equal(b$r50, d[round(N / 2)])
    expect_equal(b$ralpha, d[round(N * (1 - b$alpha))])
    expect_true(all(b$dist[b$outlier_index] > b$whisker_dis))
  }
})

test_that("`boxplot_region()` tests the whiskers last", {
  # the whisker may fall inside ralpha; beyond it is always an outlier
  reg <- boxplot_region(c(0.5, 1.5, 2.5, 4), c(1, 3, 2), 0.05)
  expect_equal(reg$level, c(1L, 2L, 4L, 4L))
  expect_equal(reg$region[2], "95% central region")
})

test_that("`align_new()` returns boundary-preserving warps", {
  a <- align_new(out_median, simu_data$f[, 1:3])
  expect_equal(dim(a$fn), c(101, 3))
  expect_equal(dim(a$qn), c(101, 3))
  expect_equal(dim(a$gam), c(101, 3))
  expect_lt(max(abs(a$gam[1, ])), 1e-12)
  expect_lt(max(abs(a$gam[101, ] - 1)), 1e-12)
  a1 <- align_new(out_median, simu_data$f[, 1])
  expect_equal(a1$gam[, 1], a$gam[, 1])
  expect_error(align_new(out_median, simu_data$f[1:50, ]))
})

test_that("`classify()` outputs are well formed", {
  fnew <- simu_data$f[, 1:3]
  for (b in list(amp_box, ph_box)) {
    c <- classify(b, fnew)
    expect_length(c$dist, 3)
    expect_true(all(is.finite(c$dist)))
    expect_true(all(c$level %in% 1:4))
    expect_true(all(c$side %in% c("Q1", "Q3")))
    expect_equal(dim(c$gam), dim(fnew))
  }
})

test_that("`classify()` places the median inside and a scaled function out", {
  expect_equal(classify(amp_box, out_median$fmean)$level, 1L)
  expect_equal(classify(amp_box, 10 * simu_data$f[, 1])$level, 4L)
})
