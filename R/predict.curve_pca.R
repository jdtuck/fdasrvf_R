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
  km <- object$karcher_mean
  if (is.null(newdata)) {
    newdata <- km$beta
  }

  dims <- dim(newdata)
  if (length(dims) == 2) {
    dims <- c(dims, 1)
    dim(newdata) <- dims
  }
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]
  mu <- km$mu
  mode <- km$mode
  scale <- km$scale

  basis <- NULL
  if (mode == "C" && scale)
    basis <- find_basis_normal(mu)

  # Align each curve to the Karcher mean and compute its shooting vector the
  # same way multivariate_karcher_mean() does for the curves it was fit on
  v <- matrix(0, L * M, N)
  for (ii in 1:N) {
    q1 <- curve_to_srvf(newdata[, , ii], scale = scale)$q
    out <- find_rotation_seed_unique(
      q1 = mu,
      q2 = q1,
      mode = mode,
      alignment = km$alignment,
      rotation = km$rotation,
      scale = scale,
      lambda = km$lambda
    )
    w <- inverse_exponential_map(out$q2best, mu, scale = scale)
    if (!is.null(basis))
      w <- project_tangent(w, mu, basis)
    v[, ii] <- c(w)
  }

  t(object$U) %*% (v - object$VM)
}
