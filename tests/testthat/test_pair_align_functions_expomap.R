library(fdasrvf)
testthat::context("Pair align functions bayes expomap")

testthat::test_that("Verify pair_align_functions_expomap() works as intended", {
  set.seed(1)
  data('simu_data')
  myf1 <- simu_data$f[,1]
  myf2 <- simu_data$f[,2]
  mytime <- simu_data$time
  myzpcn <- list(betas = c(0.5, 0.05, 0.005, 0.0001), probs = c(0.1, 0.1, 0.7, 0.1))
  out <- pair_align_functions_expomap(myf1, myf2, timet = mytime, iter = 1e4,
    alpha0 = 0.1, beta0 = 0.1, zpcn = myzpcn, extrainfo = TRUE)
  # verify function hasn't been changed
  testthat::expect_equal(mean(out$gamma$y), 0.5083758020)
  testthat::expect_equal(mean(out$g.coef), 0.0009790850)
  testthat::expect_equal(mean(out$psi$y), 0.9994264515)
  testthat::expect_equal(mean(out$sigma1), 0.1190701553)
  testthat::expect_equal(mean(out$gamma_q025), 0.5047645723)
  testthat::expect_equal(mean(out$gamma_q975), 0.5121841722)
  testthat::expect_equal(mean(out$gamma_mat), 0.5083537445)
  testthat::expect_equal(sd(out$xdist), 0.0032652812)
  testthat::expect_equal(sd(out$ydist), 0.0028818982)
  # verify functions match approximately
  testthat::expect_equal(sum(simu_data$f[,1] - out$f2_warped), -10.32011970467)
  # verify acceptance rate
  testthat::expect_equal(mean(out$accept), 0.1701170117)
  testthat::expect_equal(prod(log(table(out$betas.ind))), 1832.14013976635)
  # verify extrainfo=FALSE excludes appropriate output
  outSmall <- pair_align_functions_expomap(myf1, myf2, timet = mytime, iter = 1e4,
    alpha0 = 0.1, beta0 = 0.1, extrainfo = FALSE)
  testthat::expect_false(is.null(outSmall$f2_warped))
  testthat::expect_false(is.null(outSmall$gamma))
  testthat::expect_false(is.null(outSmall$g.coef))
  testthat::expect_false(is.null(outSmall$psi))
  testthat::expect_false(is.null(outSmall$sigma1))
  testthat::expect_true(is.null(outSmall$accept))
  testthat::expect_true(is.null(outSmall$betas.ind))
  testthat::expect_false(is.null(out$betas.ind))
  testthat::expect_true(is.null(outSmall$logl))
  testthat::expect_false(is.null(out$logl))
  testthat::expect_true(is.null(outSmall$gamma_mat))
  testthat::expect_true(is.null(outSmall$gamma_q025))
  testthat::expect_true(is.null(outSmall$gamma_q975))
  testthat::expect_true(is.null(outSmall$sigma_eff_size))
  testthat::expect_false(is.null(out$sigma_eff_size))
  testthat::expect_true(is.null(outSmall$psi_eff_size))
  testthat::expect_false(is.null(out$psi_eff_size))
  testthat::expect_true(is.null(outSmall$xdist))
  testthat::expect_false(is.null(out$xdist))
  testthat::expect_true(is.null(outSmall$ydist))
  testthat::expect_false(is.null(out$ydist))
})

testthat::test_that("Verify init.coef arguments work correctly", {
  data('simu_data')
  myf1 <- simu_data$f[1:20,1]
  myf2 <- simu_data$f[1:20,2]
  mytime <- simu_data$time[1:20]
  # must work
  out <- pair_align_functions_expomap(myf1, myf2, timet = mytime, init.coef = rep(0,6))
  testthat::expect_error(pair_align_functions_expomap(myf1, myf2, timet = mytime, init.coef = rep(0, 7)),
    'Length of init.coef must be even')
  testthat::expect_error(pair_align_functions_expomap(myf1, myf2, mytime, init.coef = rep(-1, 10)),
    'Invalid initial value of g')
})

testthat::test_that("Verify pair_align_functions_expomap() throws appropriate errors", {
  data('simu_data')
  myf1 <- simu_data$f[,1]
  myf2 <- simu_data$f[,2]
  mytime <- simu_data$time
  shortf <- myf1[1:10]
  shortt <- mytime[1:10]
  testthat::expect_error(pair_align_functions_expomap(shortf, myf2, mytime, iter = 1e3),
    'Length of f1 and f2 must be equal')
  testthat::expect_error(pair_align_functions_expomap(myf1, shortf, mytime, iter = 1e3),
    'Length of f1 and f2 must be equal')
  testthat::expect_error(pair_align_functions_expomap(myf1, myf2, shortt, iter = 1e3),
    'Length of f1 and timet must be equal')
  testthat::expect_error(pair_align_functions_expomap(shortf, shortf, mytime, iter = 1e3),
    'Length of f1 and timet must be equal')
  testthat::expect_error(pair_align_functions_expomap(shortf, myf2, shortt, iter = 1e3),
    'Length of f1 and f2 must be equal')
  testthat::expect_error(pair_align_functions_expomap(myf1, shortf, shortt, iter = 1e3),
    'Length of f1 and f2 must be equal')
  testthat:: expect_error(pair_align_functions_expomap(myf1, myf2, mytime, iter = 1e3,
    zpcn = list(betas = c(0.9, 0.7), probs = c(0.5, 0.4, 0.1))),
    'In zpcn, betas must equal length of probs')
})

testthat::test_that("Check calcY() works as intended", {
  testthat::expect_equal(sum(fdasrvf:::calcY(2, 0:2)), 0.1155056305)
  testthat::expect_equal(sum(fdasrvf:::calcY(0, 0:2)), 3.0)
})

testthat::test_that("Check cuL2norm2() works as intended", {
  testthat::expect_equal(prod(fdasrvf:::cuL2norm2(0:2, 1:3)), 22.5)
  testthat::expect_equal(prod(fdasrvf:::cuL2norm2(c(1,0,2), 1:3)), 18.75)
  testthat::expect_equal(prod(fdasrvf:::cuL2norm2(c(1,2,10), 1:3)), 545/4)
})

testthat::test_that("Check trapzCpp() works as intended", {
  testthat::expect_equal(fdasrvf:::trapzCpp(0:2, 1:3), 4)
  testthat::expect_equal(fdasrvf:::trapzCpp(c(1,0,2), 1:3), 3.5)
  testthat::expect_equal(fdasrvf:::trapzCpp(c(1,2,10), 1:3), 43/2)
})

testthat::test_that("Check order_l2norm() works as intended", {
  testthat::expect_equal(fdasrvf:::order_l2norm(0:2, 1:3), 3)
  testthat::expect_equal(fdasrvf:::order_l2norm(c(1,0,2), 1:3), sqrt(15/2))
  testthat::expect_equal(fdasrvf:::order_l2norm(c(1,2,10), 1:3), sqrt(109/2))
})

