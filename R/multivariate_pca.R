#' Curve PCA
#'
#' Calculate principal directions of a set of curves
#'
#' @param align_data fdacurve object from [multivariate_karcher_mean] of aligned data
#' @param no number of components
#' @param var_exp compute no based on value percent variance explained (example: 0.95)
#'                will override `no`
#' @param ci geodesic standard deviations (default = c(-1,0,1))
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param showplot show plots of principal directions (default = TRUE)
#' @return Returns a curve_pca object containing \item{latent}{singular values}
#' \item{U}{singular vectors}
#' \item{coef}{principal coefficients}
#' \item{pd}{principal directions}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' align_data <- multivariate_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' out <- multivariate_pca(align_data)
multivariate_pca <- function(align_data,
                             no = 3,
                             var_exp = NULL,
                             ci = c(-1, 0, 1),
                             mode = "O",
                             showplot = TRUE) {
  v = align_data$v
  mu = align_data$mu
  scale = align_data$scale
  K = multivariate_karcher_cov(align_data)

  n = nrow(mu)
  T1 = ncol(mu)

  # SVD
  out = svd(K)
  U = out$u
  s = out$d

  # Parameters
  if (!is.null(var_exp)) {
    cumm_coef <- cumsum(s) / sum(s)
    tmp = which(cumm_coef <= var_exp)
    no = tmp[length(tmp)]
  }

  U = out$u[, 1:no]
  s = out$d[1:no]

  tmp = dim(v)
  N1 = tmp[3]
  VM = apply(v, c(1, 2), mean)
  VM = c(VM)

  # express shapes as coefficients
  x = matrix(0, no, N1)
  for (ii in 1:N1) {
    tmpv = c(v[, , ii])
    x[, ii] = t(U) %*% (tmpv - VM)
  }

  N = length(ci)
  pd = array(list(), c(no, N))
  for (m in 1:no) {
    for (i in 1:N) {
      v1 = VM + ci[i] * sqrt(s[m]) * U[, m]
      tmp_scale = 1

      dim(v1) = c(n, T1)
      q2n = exponential_map(v1, mu, scale)
      if (mode == "C") {
        q2n = project_curve(q2n)
      }
      p = q_to_curve(q2n, tmp_scale)

      pd[m, i][[1]] = p

    }
  }

  curve_pca = list(
    latent = s,
    U = U,
    coef = x,
    pd = pd,
    VM = VM,
    stds = ci,
    karcher_mean = align_data
  )

  class(curve_pca) <- "curve_pca"

  if (showplot) {
    plot(curve_pca)
  }

  return(curve_pca)
}
