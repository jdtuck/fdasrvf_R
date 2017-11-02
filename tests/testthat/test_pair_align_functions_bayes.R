library(fdasrvf)
context("Pair align functions bayes Cheng")

test_that("Verify pair_align_functions_bayes() works as intended", {
  set.seed(1)
  data('simu_data')
  out <- pair_align_functions_bayes(f1 = simu_data$f[,1], f2 = simu_data$f[,2],
    timet = simu_data$time, iter = 5000, times = 5, tau = 2, powera = 1,
    showplot = FALSE, extrainfo = TRUE)
  # verify functions match approximately
  expect_equal(sum(simu_data$f[,1] - out$f2_a), -10.1946817)
  expect_equal(sum(simu_data$f[,1] - out$f2_q), -10.0712115)
  # verify warning issued for times = 2
  expect_warning(pair_align_functions_bayes(f1 = simu_data$f[,1],
    f2 = simu_data$f[,2], timet = simu_data$time, iter = 1000, times = 2,
    tau = 2, powera = 1, showplot = FALSE, extrainfo = FALSE))
  # verify acceptance rate
  expect_equal(mean(out$pct_accept), 0.0304315789)
})

