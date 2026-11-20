test_that("`dpmeans_align()` learns the number of clusters", {
  t <- seq(0, 1, length.out = 101)
  f <- cbind(
    sapply(c(0.35, 0.4, 0.45), function(m) dnorm(t, m, 0.07)),
    sapply(c(0, 0.1, 0.2), function(s) dnorm(t, 0.3, 0.05) +
             dnorm(t, 0.7, 0.05) + s)
  )
  out <- dpmeans_align(f, t, lambda = 2)
  expect_s3_class(out, "fdakma")
  expect_equal(out$labels, c(1, 1, 1, 2, 2, 2))
  expect_equal(dim(out$templates), c(1, 101, 2))
  expect_length(out$fn, 2)
  expect_length(out$gamI, 2)
  expect_equal(dim(out$gam[[2]]), c(101, 3))

  # a huge penalty keeps a single cluster
  out1 <- dpmeans_align(f, t, lambda = 1e6)
  expect_equal(out1$labels, rep(1, 6))
  expect_equal(out1$qun[1], sum(out1$distances_to_center^2) + 1e12)
})
