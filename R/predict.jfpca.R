#' Elastic Prediction for functional PCA
#'
#' This function performs projection of new functions on fPCA basis
#'
#' @param object Object of class inheriting from "jointFPCA"
#' @param newdata An optional matrix in which to look for functions with which to predict. If omitted, the original functions are used.
#' @param ... additional arguments affecting the predictions produced
#' @return Returns a matrix
#' \item{a}{principle coefficients}
#' @keywords srvf alignment regression
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
predict.jfpca <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    newdata = object$warp_data$f0
  }
  q1 = f_to_srvf(newdata, object$warp_data$time)
  M = length(object$warp_data$time)
  N = ncol(newdata)
  gam = matrix(0, M, N)
  fn = matrix(0, M, N)
  qn = matrix(0, M, N)
  for (ii in 1:N) {
    gam[, ii] = optimum.reparam(
      object$warp_data$mqn,
      object$warp_data$time,
      q1[, ii],
      object$warp_data$time,
      method = object$warp_data$call$optim_method
    )
    fn[, ii] = warp_f_gamma(newdata[, ii], object$warp_data$time, gam[, ii])
    qn[, ii] = f_to_srvf(fn[, ii], object$warp_data$time)
  }

  m_new <- sign(fn[object$id, ]) * sqrt(abs(fn[object$id, ]))  # scaled version
  qn1 <- rbind(qn, m_new)

  no = ncol(object$U)
  psi = matrix(0, M, N)
  vec = matrix(0, M, N)
  time = np.linspace(0, 1, M)
  binsize <- mean(diff(time))
  for (i in 1:N) {
    psi[, i] = sqrt(gradient(gam[, i], binsize))
    vec[, i] <- inv_exp_map(object$mu_psi, psi[, i])
  }

  g <- rbind(qn1, object$C * vec)
  a <- matrix(0, N, no)
  for (i in 1:N) {
    for (j in 1:no) {
      a[i, j] <- (g[, i] - object$mu_g) %*% object$U[, j]
    }
  }

  return(a)
}
