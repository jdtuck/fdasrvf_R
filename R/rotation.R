# return the lower-triangular indices of a m x m matrix
LTindgen <- function(m) {
  aa <- unlist(lapply(1:(m - 1), function(x)
    (x + 1):m))
  bb <- rep(1:(m - 1), (m - 1):1)
  unname(cbind(aa, bb))
}

# principal matrix logarithm
matlog <- function(X,
                   eig = eigen(X),
                   skewsymmetric = F,
                   Retol = 1e-8) {
  xx <- eig$vectors %*% diag(log(eig$values)) %*% Conj(t(eig$vectors))
  if (max(abs(Im(xx))) < Retol)
    xx <- Re(xx)
  else{
    print("complex logarithm")
  }
  if (skewsymmetric)
    return((xx - t(xx)) / 2)
  return(xx)
}

# unique matrix exponential
matexp <- function(X, eig = eigen(X), Retol = 1e-8) {
  xx <- Re(eig$vectors %*% diag(exp(eig$values)) %*% Conj(t(eig$vectors)))
  if (max(abs(Im(xx))) < Retol)
    Re(xx)
  else{
    print("complex exponential")
    xx
  }
}

# karcher mean of rotation matrices
kmean_NDrotations <- function(Os,
                              M = Os[[1]],
                              tau = 1,
                              maxit = 100,
                              tol = 1e-8) {
  for (k in 1:maxit) {
    As <- lapply(Os, function(O)
      matlog(O %*% t(M), skewsymmetric = T))
    A <- Reduce("+", As) / length(Os)
    M <- matexp(tau * A) %*% M
    if (sum(abs(A)) < tol)
      break
  }
  info <- list(
    nit = k,
    convergence = (k != maxit),
    crit = sum(abs(A))
  )
  list(kmean = M, As = As, info = info)
}

#' Rotation Principal Component Analysis
#'
#' This function calculates functional principal component analysis
#' on rotation data from
#'
#' @param align_data fdacurve object from [multivariate_karcher_mean] of aligned data
#' @param no number of principal components to extract
#' @param var_exp compute no based on value percent variance explained (example: 0.95)
#'                will override `no`
#' @return Returns a rotpca object containing
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' @export
rotation_pca <- function(align_data,
                         no = 3,
                         var_exp = NULL) {
  N <- dim(align_data$R)[3]
  Nd <- dim(align_data$beta)[1]
  Os <- lapply(1:N, function(x)
    align_data$R[, , x]) #this transforms array to list
  kmeanOs <- kmean_NDrotations(Os) #find Karcher mean of rotations
  ltind <- LTindgen(Nd)
  Nr <- nrow(ltind)
  Rs <- sapply(kmeanOs$As, function(x)
    x[ltind])

  if (ndims(Rs) == 0) {
    Rs = t(t(Rs))
  }

  K = stats::cov(Rs)

  out = svd(K)
  s = out$d
  stdS = sqrt(s)
  U = out$u

  if (!is.null(var_exp)) {
    cumm_coef <- cumsum(s) / sum(s)
    tmp = which(cumm_coef <= var_exp)
    no = tmp[length(tmp)]
  }

  if (no > length(s)) {
    no = length(s)
  }


  Rm = colMeans(Rs)
  c = matrix(0, N, no)
  for (k in no) {
    for (i in 1:N) {
      c[i, k] = sum((Rs[i] - Rm) * U[, k])
    }
  }

  rfpca <- list()
  rfpca$latent <- s[1:no]
  rfpca$coef <- c
  rfpca$U <- U[, 1:no]
  rfpca$mean <- Rm
  rfpca$align_data <- align_data

  class(rfpca) <- "rotpca"

  return(rfpca)

}


# 1.3. extract rotation from euclidean
# Mr,Xr: Karcher mean and vector-space mapping of rotation matrices
# rotmatR <- apply(Xr,2,function(xx){
#   rr <- matrix(0,Nd,Nd)
#   rr[ltind] <- xx
#   matexp(rr - t(rr))%*%Mr #from tangent map to rotation manifold
# })
# dim(rotmatR) <- c(Nd,Nd,N)
