if (requireNamespace("testthat", quietly = TRUE)){
  library(fdasrvf)

  testthat::test_check("fdasrvf")
  testthat::snapshot_accept(files = NULL, path = "testthat")
} else{
  stop("Need testthat to run checks")
}

