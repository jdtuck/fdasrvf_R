#' Elastic Prediction for curve PCA
#'
#' This function performs projection of new curves on fPCA basis
#'
#' @param object Object of class inheriting from "curve_pca"
#' @param newdata An optional matrix in which to look for functions with which to predict. If omitted, the original functions are used.
#' @param ... additional arguments affecting the predictions produced
#' @return Returns a matrix
#' \item{a}{principal coefficients}
#' @keywords srvf alignment regression
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
predict.curve_pca <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    newdata = object$karcher_mean$beta
  }

  N = dim(newdata)[3]
  M = dim(newdata)[1]*dim(newdata)[2]
  mu = object$karcher_mean$mu
  rotated = object$karcher_mean$rotated
  mode = object$karcher_mean$mode
  lambda = object$karcher_mean$lambda
  ms = object$karcher_mean$ms
  if (mode == "C")
    basis <- find_basis_normal(mu)
  v = matrix(0, M, N)
  for (ii in 1:N) {
    q1 = curve_to_q(newdata[,,ii], object$karcher_mean$scale)$q
    out = karcher_calc(q1, mu, basis, rotated, mode, lambda, ms)
    v[, ii] = c(out$v)
  }

  no = ncol(object$U)

  a <- matrix(0, no, N)
  for (i in 1:no) {
    for (j in 1:N) {
      a[i, j] <- (v[, j] - object$VM) %*% object$U[, i]
    }
  }

  return(a)
}
