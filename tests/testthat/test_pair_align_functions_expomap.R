library(fdasrvf)
context("Pair align functions bayes expomap")

test_that("Verify pair_align_functions_expomap() works as intended", {
  set.seed(1)
  data('simu_data')
  myf1 <- simu_data$f[,1]
  myf2 <- simu_data$f[,2]
  mytime <- simu_data$time
  myzpcn <- list(betas = c(0.5, 0.05, 0.005, 0.0001), probs = c(0.1, 0.1, 0.7, 0.1))
  out <- pair_align_functions_expomap(myf1, myf2, timet = mytime, iter = 1e4,
    alpha0 = 0.1, beta0 = 0.1, zpcn = myzpcn, extrainfo = TRUE)
  # verify function hasn't been changed
  expect_equal(mean(out$gamma$y), 0.5083758020)
  expect_equal(mean(out$g.coef), 0.0009790850)
  expect_equal(mean(out$psi$y), 0.9994264515)
  expect_equal(mean(out$sigma1), 0.1190701553)
  expect_equal(mean(out$gamma_q025), 0.5047645723)
  expect_equal(mean(out$gamma_q975), 0.5121841722)
  expect_equal(mean(out$gamma_mat), 0.5083537445)
  # verify functions match approximately
  expect_equal(sum(simu_data$f[,1] - out$f2_warped), -10.32011970467)
  # verify acceptance rate
  expect_equal(mean(out$accept), 0.1701170117)
  expect_equal(prod(log(table(out$betas.ind))), 1832.14013976635)
  # verify extrainfo=FALSE excludes appropriate output
  outSmall <- pair_align_functions_expomap(myf1, myf2, timet = mytime, iter = 1e4,
    alpha0 = 0.1, beta0 = 0.1, extrainfo = FALSE)
  expect_false(is.null(outSmall$f2_warped))
  expect_false(is.null(outSmall$gamma))
  expect_false(is.null(outSmall$g.coef))
  expect_false(is.null(outSmall$psi))
  expect_false(is.null(outSmall$sigma1))
  expect_true(is.null(outSmall$accept))
  expect_true(is.null(outSmall$betas.ind))
  expect_false(is.null(out$betas.ind))
  expect_true(is.null(outSmall$logl))
  expect_false(is.null(out$logl))
  expect_true(is.null(outSmall$gamma_mat))
  expect_true(is.null(outSmall$gamma_q025))
  expect_true(is.null(outSmall$gamma_q975))
  expect_true(is.null(outSmall$sigma_eff_size))
  expect_false(is.null(out$sigma_eff_size))
  expect_true(is.null(outSmall$psi_eff_size))
  expect_false(is.null(out$psi_eff_size))
})

test_that("Verify init.coef arguments work correctly", {
  data('simu_data')
  myf1 <- simu_data$f[1:20,1]
  myf2 <- simu_data$f[1:20,2]
  mytime <- simu_data$time[1:20]
  # must work
  out <- pair_align_functions_expomap(myf1, myf2, timet = mytime, init.coef = rep(0,6))
  expect_error(pair_align_functions_expomap(myf1, myf2, timet = mytime, init.coef = rep(0, 7)), 
    'Length of init.coef must be even')
  expect_error(pair_align_functions_expomap(myf1, myf2, mytime, init.coef = rep(-1, 10)),
    'Invalid initial value of g')
})

test_that("Verify pair_align_functions_expomap() throws appropriate errors", {
  data('simu_data')
  myf1 <- simu_data$f[,1]
  myf2 <- simu_data$f[,2]
  mytime <- simu_data$time
  shortf <- myf1[1:10]
  shortt <- mytime[1:10]
  expect_error(pair_align_functions_expomap(shortf, myf2, mytime, iter = 1e3),
    'Length of f1 and f2 must be equal')
  expect_error(pair_align_functions_expomap(myf1, shortf, mytime, iter = 1e3),
    'Length of f1 and f2 must be equal')
  expect_error(pair_align_functions_expomap(myf1, myf2, shortt, iter = 1e3),
    'Length of f1 and timet must be equal')
  expect_error(pair_align_functions_expomap(shortf, shortf, mytime, iter = 1e3),
    'Length of f1 and timet must be equal')
  expect_error(pair_align_functions_expomap(shortf, myf2, shortt, iter = 1e3),
    'Length of f1 and f2 must be equal')
  expect_error(pair_align_functions_expomap(myf1, shortf, shortt, iter = 1e3),
    'Length of f1 and f2 must be equal')
  expect_error(pair_align_functions_expomap(myf1, myf2, mytime, iter = 1e3,
    zpcn = list(betas = c(0.9, 0.7), probs = c(0.5, 0.4, 0.1))),
    'In zpcn, betas must equal length of probs')
})

test_that("Check calcY() works as intended", {
  expect_equal(sum(fdasrvf:::calcY(2, 0:2)), 0.1155056305)
  expect_equal(sum(fdasrvf:::calcY(0, 0:2)), 3.0)
})

test_that("Check cuL2norm2() works as intended", {
  expect_equal(prod(fdasrvf:::cuL2norm2(0:2, 1:3)), 22.5)
  expect_equal(prod(fdasrvf:::cuL2norm2(c(1,0,2), 1:3)), 18.75)
  expect_equal(prod(fdasrvf:::cuL2norm2(c(1,2,10), 1:3)), 545/4)
})

test_that("Check trapzCpp() works as intended", {
  expect_equal(fdasrvf:::trapzCpp(0:2, 1:3), 4)
  expect_equal(fdasrvf:::trapzCpp(c(1,0,2), 1:3), 3.5)
  expect_equal(fdasrvf:::trapzCpp(c(1,2,10), 1:3), 43/2)
})

test_that("Check order_l2norm() works as intended", {
  expect_equal(fdasrvf:::order_l2norm(0:2, 1:3), 3)
  expect_equal(fdasrvf:::order_l2norm(c(1,0,2), 1:3), sqrt(15/2))
  expect_equal(fdasrvf:::order_l2norm(c(1,2,10), 1:3), sqrt(109/2))
})

