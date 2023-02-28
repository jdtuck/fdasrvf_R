test_that("`boxplot()` works", {
  out <- time_warping(simu_data$f, simu_data$time)
  amp_data <- boxplot(out, variability_type = "amplitude", what = "stats")
  ph_data <- boxplot(out, variability_type = "phase", what = "stats")

  # Amplitude boxplot
  expect_equal(class(amp_data), "ampbox")
  expect_equal(length(amp_data), 22)
  expect_equal(names(amp_data), c("median_y", "Q1", "Q3", "Q1a", "Q3a", "minn",
                                 "maxx", "outlier_index", "fmedian", "time",
                                 "qmedian", "Q1_q", "Q1a_q", "Q3_q", "Q3a_q",
                                 "min_q", "max_q", "dist", "Q1a_index",
                                 "Q1_index", "Q3a_index", "Q3_index"))
  expect_equal(length(amp_data$median_y), 101)
  expect_equal(length(amp_data$Q1), 101)
  expect_equal(length(amp_data$Q3), 101)
  expect_equal(length(amp_data$Q1a), 101)
  expect_equal(length(amp_data$Q3a), 101)
  expect_equal(length(amp_data$minn), 101)
  expect_equal(length(amp_data$maxx), 101)
  expect_true(is.null(amp_data$outlier_index))
  expect_equal(length(amp_data$fmedian), 101)
  expect_equal(length(amp_data$time), 101)
  expect_equal(length(amp_data$qmedian), 101)
  expect_equal(length(amp_data$Q1_q), 101)
  expect_equal(length(amp_data$Q1a_q), 101)
  expect_equal(length(amp_data$Q3_q), 101)
  expect_equal(length(amp_data$Q3a_q), 101)
  expect_equal(length(amp_data$min_q), 101)
  expect_equal(length(amp_data$max_q), 101)
  expect_equal(length(amp_data$dist), 21)
  expect_equal(amp_data$Q1a_index, 3)
  expect_equal(amp_data$Q1_index, 3)
  expect_equal(amp_data$Q3a_index, 21)
  expect_equal(amp_data$Q3_index, 5)
  expect_snapshot(amp_data)
})
