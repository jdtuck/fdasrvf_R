#' Long Run Covariance Matrix Estimation for Multivariate Time Series
#'
#' This function estimates the long run covariance matrix of a given multivariate data sample.
#'
#' @param mdobj A multivariate data object
#' @param h The bandwidth parameter. It is strictly non-zero. Choosing the bandwidth parameter to be zero is identical
#' to estimating covariance matrix assuming iid data.
#' @param kern_type Kernel function to be used for the estimation of the long run covariance
#' matrix. The choices are \code{c("BT", "PR", "SP", "FT")} which are respectively, bartlett, parzen, simple and flat-top kernels.
#' By default the function uses a \code{"barlett"} kernel.
#' @return Returns long run covariance matrix

# this is for the computation of Long Run Variance of \Theta
LongRunCovMatrix <- function(mdobj, h=0, kern_type = "bartlett"){
  N = ncol(mdobj)
  D = nrow(mdobj)
  Kernel <- function(i, h) {
    x = i/h
    if (kern_type == "flat") {
      return(1)
    }
    if (kern_type == "simple") {
      return(0)
    }
    if (kern_type == "bartlett") {
      return(1 - x)
    }
    if (kern_type == "flat_top") {
      if (x < 0.1) {
        return(1)
      } else {
        if (x >= 0.1 & x < 1.1) {
          return(1.1 - x)
        } else {
          return(0)
        }
      }
    }
    if (kern_type == "parzen") {
      if (x < 1/2) {
        return(1 - 6 * x^2 + 6 * abs(x)^3)
      } else {
        return(2 * (1 - abs(x))^3)
      }
    }
  }
  D_mat = matrix(0, D, D)
  cdata = mdobj
  # Long Run Cov Est
  for (k in 1:D) {
    for (r in k:D) {
      s = cdata[k, 1:N] %*% cdata[r, 1:N]
      if (h > 0) {
        for (i in 1:h) {
          a = cdata[k, 1:(N - i)] %*% cdata[r,(i + 1):N]
          a = a + cdata[r, 1:(N - i)] %*% cdata[k,(i + 1):N]
          s = s + Kernel(i, h) * a
        }
      }
      D_mat[k, r] = s
      D_mat[r, k] = D_mat[k, r]
    }
  }
  D_mat/N
}
