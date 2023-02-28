test_that("`align_fPCA()` works", {
  out <- align_fPCA(simu_data$f, simu_data$time, showplot = FALSE)
  M <- dim(simu_data$f)[1]
  N <- dim(simu_data$f)[2]
  n_sd <- 2
  n_components <- 3
  expect_equal(dim(out$f0), c(M, N))
  expect_equal(dim(out$q0), c(M, N))
  expect_equal(dim(out$fn), c(M, N))
  expect_equal(dim(out$qn), c(M, N))
  expect_equal(length(out$mqn), M)
  expect_equal(dim(out$gam), c(M, N))
  expect_equal(length(out$vfpca), 5)
  expect_equal(names(out$vfpca), c("q_pca", "f_pca", "latent", "coef", "U"))
  expect_equal(dim(out$vfpca$q_pca), c(M + 1, 2 * n_sd + 1, n_components))
  expect_equal(dim(out$vfpca$f_pca), c(M, 2 * n_sd + 1, n_components))
  expect_equal(length(out$vfpca$latent), M + 1)
  expect_equal(dim(out$vfpca$coef), c(N, n_components))
  expect_equal(dim(out$vfpca$U), c(M + 1, n_components))
  expect_equal(length(out$Dx), 51)
  expect_snapshot(out)
})
