# Smoke / return-structure tests for the heavier modeling routines that had no
# direct test. These functions are exercised on small, seeded inputs with
# plotting disabled where possible; the assertions check the documented return
# structure and key dimensions rather than exact numeric values (mirroring the
# style of test-kmeans.R and test-align-fpca.R).
#
# Note: the elastic (m)(l)pcr regression fitters are already exercised by
# test-predict-pcr.R, and gauss_model() by test-regressions.R, so they are not
# duplicated here.

test_that("`multivariate_karcher_cov()` returns a symmetric covariance matrix", {
  km <- multivariate_karcher_mean(fdasrvf::beta[, , 1, 1:4], maxit = 2)
  K <- multivariate_karcher_cov(km)
  expect_true(is.matrix(K))
  expect_equal(nrow(K), ncol(K))
  # A covariance matrix is symmetric.
  expect_equal(K, t(K), tolerance = 1e-8)
})

test_that("`rotation_pca()` returns a rotpca object", {
  km <- multivariate_karcher_mean(
    fdasrvf::beta[, , 1, 1:5], rotation = TRUE, scale = TRUE, maxit = 2
  )
  out <- rotation_pca(km, no = 2)
  expect_s3_class(out, "rotpca")
  expect_true(all(c("latent", "coef", "U", "mean") %in% names(out)))
})

test_that("`outlier.detection()` returns the flagged outlier SRVFs", {
  q_outlier <- outlier.detection(
    q = toy_warp$q0,
    time = toy_data$time,
    mq = toy_warp$mqn,
    k = 0.1
  )
  # The function returns the columns of `q` flagged as outliers (subset of the
  # grid-by-curve matrix); `k = 0.1` is a loose threshold so at least one is
  # returned. Coerce to a matrix since a single outlier drops to a vector.
  q_outlier <- as.matrix(q_outlier)
  expect_equal(nrow(q_outlier), nrow(toy_warp$q0))
  expect_true(ncol(q_outlier) <= ncol(toy_warp$q0))
  expect_true(all(is.finite(q_outlier)))
})

test_that("`joint_gauss_model()` generates M x n aligned samples", {
  withr::with_seed(1234, {
    out <- joint_gauss_model(simu_warp, n = 3)
  })
  M <- length(simu_warp$time)
  expect_equal(dim(out$fs), c(M, 3))
  expect_false(anyNA(out$fs))
})

test_that("`elastic_amp_change_ff()` detects an amplitude change point", {
  withr::with_seed(1234, {
    out <- elastic_amp_change_ff(
      simu_data$f, simu_data$time, d = 50, showplot = FALSE
    )
  })
  expect_true(all(c("pvalue", "change") %in% names(out)))
  expect_true(out$pvalue >= 0 && out$pvalue <= 1)
  expect_true(out$change >= 1 && out$change <= ncol(simu_data$f))
})

test_that("`elastic_ph_change_ff()` detects a phase change point", {
  withr::with_seed(1234, {
    out <- elastic_ph_change_ff(
      simu_data$f, simu_data$time, d = 50, showplot = FALSE
    )
  })
  expect_true(all(c("pvalue", "change") %in% names(out)))
  expect_true(out$pvalue >= 0 && out$pvalue <= 1)
})

test_that("`elastic_change_fpca()` detects a change point via fPCA", {
  withr::with_seed(1234, {
    out <- elastic_change_fpca(
      simu_data$f, simu_data$time, d = 50, showplot = FALSE
    )
  })
  expect_true(all(c("pvalue", "change") %in% names(out)))
  expect_true(out$pvalue >= 0 && out$pvalue <= 1)
})
