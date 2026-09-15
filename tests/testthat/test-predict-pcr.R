f <- simu_data$f
time <- simu_data$time
newdata <- f[, 1:5]

# Records the arguments each predict method hands to `optimum.reparam()`
capture_reparam_args <- function(code) {
  orig <- optimum.reparam
  calls <- list()
  local_mocked_bindings(optimum.reparam = function(...) {
    args <- list(...)
    calls[[length(calls) + 1]] <<- args[c("lambda", "pen", "method")]
    orig(...)
  })
  force(code)
  calls
}

# elastic.lpcr.regression() and elastic.mlpcr.regression() always align with
# `parallel = TRUE`, which spawns more processes than R CMD check allows, so
# align serially instead
local_serial_time_warping <- function(env = parent.frame()) {
  orig <- time_warping
  local_mocked_bindings(
    time_warping = function(f, time, ..., parallel) orig(f, time, ...),
    .env = env
  )
}

set.seed(1)
y_pcr <- colMeans(f) + stats::rnorm(ncol(f), sd = 0.1)
fit_pcr <- suppressMessages(
  elastic.pcr.regression(f, y_pcr, time, pca.method = "vert", no = 3)
)

test_that("`predict.pcr()` predicts newdata", {
  out <- predict(fit_pcr, newdata = newdata, y = y_pcr[1:5])
  expect_length(out$y_pred, 5)
  expect_true(all(is.finite(out$y_pred)))
  expect_equal(out$SSE, sum((y_pcr[1:5] - out$y_pred)^2))
})

test_that("`predict.pcr()` predicts newdata for every pca method", {
  for (pca.method in c("combined", "horiz")) {
    fit <- suppressMessages(
      elastic.pcr.regression(f, y_pcr, time, pca.method = pca.method, no = 3)
    )
    out <- predict(fit, newdata = newdata)
    expect_length(out$y_pred, 5)
    expect_true(all(is.finite(out$y_pred)))
  }
})

test_that("`predict.pcr()` aligns newdata with the settings used to fit", {
  fit <- fit_pcr
  fit$warp_data$call$lambda <- 0.01
  fit$warp_data$call$penalty_method <- "l2gam"
  fit$warp_data$call$optim_method <- "RBFGS"
  calls <- capture_reparam_args(predict(fit, newdata = newdata))
  expect_length(calls, ncol(newdata))
  for (args in calls)
    expect_equal(args, list(lambda = 0.01, pen = "l2gam", method = "RBFGS"))
})

test_that("`predict.lpcr()` predicts newdata", {
  local_serial_time_warping()
  y <- rep(c(-1, 1), length.out = ncol(f))
  fit <- suppressMessages(
    elastic.lpcr.regression(f, y, time, pca.method = "vert", no = 3)
  )
  calls <- capture_reparam_args(
    out <- predict(fit, newdata = newdata, y = y[1:5])
  )
  expect_length(out$y_pred, 5)
  expect_true(all(out$y_pred >= 0 & out$y_pred <= 1))
  expect_true(all(out$y_labels %in% c(-1, 1)))
  expect_true(out$PC >= 0 && out$PC <= 1)
  expect_equal(calls[[1]], list(lambda = 0, pen = "roughness", method = "DP"))
})

test_that("`predict.mlpcr()` predicts newdata", {
  local_serial_time_warping()
  y <- rep(1:3, length.out = ncol(f))
  fit <- suppressMessages(
    elastic.mlpcr.regression(f, y, time, pca.method = "vert", no = 3)
  )
  calls <- capture_reparam_args(
    out <- predict(fit, newdata = newdata, y = y[1:5])
  )
  expect_equal(dim(out$y_pred), c(5, 3))
  expect_true(all(out$y_labels %in% 1:3))
  expect_length(out$PC, 3)
  expect_true(out$PC.comb >= 0 && out$PC.comb <= 1)
  expect_equal(calls[[1]], list(lambda = 0, pen = "roughness", method = "DP"))
})

test_that("fpca predict methods align newdata with the settings used to fit", {
  warp <- suppressMessages(time_warping(f, time, max_iter = 1))
  warp$call$lambda <- 0.01
  warp$call$penalty_method <- "l2gam"
  warp$call$optim_method <- "RBFGS"
  expected <- list(lambda = 0.01, pen = "l2gam", method = "RBFGS")
  fits <- list(
    vfpca = vertFPCA(warp, no = 3, showplot = FALSE),
    hfpca = horizFPCA(warp, no = 3, showplot = FALSE),
    jfpca = jointFPCA(warp, no = 3, showplot = FALSE),
    jfpcah = jointFPCAh(warp, showplot = FALSE)
  )
  for (nm in names(fits)) {
    calls <- capture_reparam_args(out <- predict(fits[[nm]], newdata))
    expect_length(calls, ncol(newdata))
    expect_equal(calls[[1]], expected, label = nm)
    expect_equal(nrow(out), ncol(newdata), label = nm)
  }
})

test_that("`predict.jfpcah()` reproduces the fitted coefficients on training data", {
  warp <- suppressMessages(time_warping(f, time, max_iter = 1))
  fit <- jointFPCAh(warp, showplot = FALSE)
  pred <- predict(fit)
  expect_equal(dim(pred), dim(fit$coef))
  # time_warping() and predict() discretize the aligned SRSFs differently, so
  # the match is close rather than exact
  expect_lt(norm(pred - fit$coef, "F") / norm(fit$coef, "F"), 0.12)
})
