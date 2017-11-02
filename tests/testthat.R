if (requireNamespace("testthat", quietly = TRUE)){
  library(fdasrvf)

  testthat::test_check("fdasrvf")
} else{
  stop("Need testthat to run checks")
}

