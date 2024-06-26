#' Elastic Amplitude Changepoint Detection
#'
#' This function identifies a amplitude changepoint using a fully functional
#' approach
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param d number of monte carlo iterations of Brownian Bridge (default = 1000)
#' @param h window selection of long range covariance function (default = 0)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param showplot show results plots (default = T)
#' @return Returns a list object containing
#' \item{pvalue}{p value}
#' \item{change}{indice of changepoint}
#' \item{DataBefore}{functions before changepoint}
#' \item{DataAfter}{functions after changepoint}
#' \item{MeanBefore}{mean function before changepoint}
#' \item{MeanAfter}{mean function after changepoint}
#' \item{WarpingBefore}{warping functions before changepoint}
#' \item{WarpingAfter}{warping functions after changepoint}
#' \item{WarpingMeanBefore}{mean warping function before changepoint}
#' \item{WarpingMeanAfter}{mean warping function after changepoint}
#' \item{change_fun}{amplitude change function}
#' \item{Sn}{test statistic values}
#' \item{mu}{mean srsfs}
#' \item{mu_f}{mean functions}
#' @keywords srvf alignment changepoint
#' @references J. D. Tucker and D. Yarger, “Elastic Functional Changepoint
#'   Detection of Climate Impacts from Localized Sources”, Envirometrics,
#'   10.1002/env.2826, 2023.
#' @export
elastic_amp_change_ff <- function(f, time, d = 1000, h = 0, smooth_data=FALSE, sparam=25, showplot = TRUE) {

  if (smooth_data) f <- smooth.data(f, sparam)

  out <- time_warping(f, time, parallel = T)

  N1 <- dim(f)[2]
  M <- dim(f)[1]
  mu <- matrix(0, nrow = M, ncol = N1)
  mu_f <- matrix(0, nrow = M, ncol = N1)

  # Compute Karcher mean for first i+1 functions
  mu[, 1] <- out$qn[,1]
  mu_f[, 1] <- out$fn[,1]
  for (i in 2:N1) {
    mu[, i] <- rowSums(out$qn[,1:i])
    mu_f[, i] <- rowSums(out$fn[,1:i])
  }

  # compute test statistic
  Sn <- (1:N1)
  Sn[1] <- 0
  for (j in (2:N1)) {
    Sn[j] <- 1/M * pvecnorm(1/sqrt(N1)*(mu_f[,j]-(j/N1)*mu_f[,N1]))^2
  }
  k.star <- min(which(Sn == max(Sn)))
  Tn <- max(Sn)

  # compute means on either side of the changepoint
  dat.a <- f[, 1:k.star]
  dat.b <- f[, (k.star + 1):N1]
  warp.a <- out$warping_functions[,1:k.star]
  warp.b <- out$warping_functions[,(k.star + 1):N1]
  mean.a <- rowMeans(out$fn[,1:k.star])
  mean.b <- rowMeans(out$fn[,(k.star + 1):N1])
  out.gam.mean <- SqrtMean(warp.a)
  warp_mean_a = out.gam.mean$gam_mu
  out.gam.mean <- SqrtMean(warp.b)
  warp_mean_b = out.gam.mean$gam_mu
  delta <- mean.a - mean.b

  # center your data
  centered_data = matrix(0, M, N1)
  for (i in (1:(k.star-1))){
    centered_data[,i] = out$fn[,i] - mean.a
  }
  for (i in (k.star:N1)){
    centered_data[,i] = out$fn[,i] - mean.b
  }

  D_mat <- LongRunCovMatrix(centered_data, h = h)
  eigen_struct <- eigen(D_mat, symmetric = TRUE)
  lambda <- Re(eigen_struct$values)/M

  asymp <- function(N) {
    BridgeLam <- matrix(0, M, N)
    for (j in (1:M)) {
      BridgeLam[j, ] <- lambda[j] * (BBridge(0, 0, 0, 1, N - 1)^2)
    }
    max(colSums(BridgeLam))
  }

  Values <- sapply(1:d, function(k) asymp(N1))
  z <- Tn <= Values
  p <- length(z[z == TRUE]) / length(z)

  # Plot
  if (showplot == TRUE) {
    graphics::par(mfrow = c(1, 3))
    graphics::matplot(f, type = "l", col = "grey", main = "Functional Data", ylab = "values")
    for (i in 1:ncol(dat.a)) {
      graphics::lines(dat.a[, i], col = "pink")
    }
    for (i in 1:ncol(dat.b)) {
      graphics::lines(dat.b[, i], col = "lightblue")
    }
    graphics::lines(mean.b, col = "blue")
    graphics::lines(mean.a, col = "red")
    graphics::legend("topleft", c("before", "after"), col = c("blue", "red"), lty = c(1, 1), cex = 0.5)
    graphics::plot(delta, type = "l", main = "Estimated Change Function", ylab = "values")

    graphics::matplot(out$warping_functions, type = "l", col = "grey", main = "Warping Functions", ylab = "values")
    for (i in 1:ncol(warp.b)) {
      graphics::lines(warp.b[, i], col = "pink")
    }
    for (i in 1:ncol(warp.a)) {
      graphics::lines(warp.a[, i], col = "lightblue")
    }
    graphics::lines(warp_mean_a, col = "blue")
    graphics::lines(warp_mean_b, col = "red")
    graphics::legend("topleft", c("before", "after"), col = c("blue", "red"), lty = c(1, 1), cex = 0.5)
  }

  out <- list(
    pvalue = p, change = k.star,
    DataBefore = dat.a, DataAfter = dat.b,
    MeanBefore = mean.a, MeanAfter = mean.b,
    WarpingMeanBefore = warp_mean_a, WarpingMeanAfter = warp_mean_b,
    WarpingBefore = warp.a, WarpingAfter = warp.b,
    change_fun = delta, Sn = Sn, mu=mu, mu_f=mu_f
  )

  return(out)
}


#' Elastic Phase Changepoint Detection
#'
#' This function identifies a phase changepoint using a fully functional
#' approach
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param d number of monte carlo iterations of Brownian Bridge (default = 1000)
#' @param h window selection of long range covariance function (default = 0)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param showplot show results plots (default = T)
#' @return Returns a list object containing
#' \item{pvalue}{p value}
#' \item{change}{indice of changepoint}
#' \item{DataBefore}{functions before changepoint}
#' \item{DataAfter}{functions after changepoint}
#' \item{MeanBefore}{mean function before changepoint}
#' \item{MeanAfter}{mean function after changepoint}
#' \item{WarpingBefore}{warping functions before changepoint}
#' \item{WarpingAfter}{warping functions after changepoint}
#' \item{WarpingMeanBefore}{mean warping function before changepoint}
#' \item{WarpingMeanAfter}{mean warping function after changepoint}
#' \item{change_fun}{amplitude change function}
#' \item{Sn}{test statistic values}
#' \item{mu}{mean shooting vectors}
#' @keywords srvf alignment changepoint
#' @references J. D. Tucker and D. Yarger, “Elastic Functional Changepoint
#'   Detection of Climate Impacts from Localized Sources”, Envirometrics,
#'   10.1002/env.2826, 2023.
#' @export
elastic_ph_change_ff <- function(f, time, d = 1000, h = 0, smooth_data=FALSE, sparam=25, showplot = TRUE) {

  if (smooth_data) f <- smooth.data(f, sparam)

  out <- time_warping(f, time, parallel = T)

  N1 <- dim(f)[2]
  M <- dim(f)[1]
  mu <- matrix(0, nrow = M, ncol = N1)

  # Compute Karcher mean of warping functions
  out.mean = SqrtMean(out$warping_functions)
  vec = out.mean$vec

  mu[, 1] <- vec[,1]
  for (i in 2:N1) {
    mu[, i] <- rowSums(vec[,1:i])
  }

  # compute test statistic
  Sn <- (1:N1)
  Sn[1] <- 0
  for (j in (2:N1)) {
    Sn[j] <- 1/M * pvecnorm(1/sqrt(N1)*(mu[,j]-(j/N1)*mu[,N1]))^2
  }
  k.star <- min(which(Sn == max(Sn)))
  Tn <- max(Sn)

  # compute means on either side of the changepoint
  dat.a <- f[, 1:k.star]
  dat.b <- f[, (k.star + 1):N1]
  warp.a <- out$warping_functions[,1:k.star]
  warp.b <- out$warping_functions[,(k.star + 1):N1]
  mean.a <- rowMeans(out$fn[,1:k.star])
  mean.b <- rowMeans(out$fn[,(k.star + 1):N1])
  out.gam.mean <- SqrtMean(warp.a)
  warp_mean_a = out.gam.mean$gam_mu
  mu.a = rowMeans(out.gam.mean$vec)
  out.gam.mean <- SqrtMean(warp.b)
  warp_mean_b = out.gam.mean$gam_mu
  mu.b = rowMeans(out.gam.mean$vec)
  delta <- mean.a - mean.b

  # center your data
  centered_data = matrix(0, M, N1)
  for (i in (1:(k.star-1))){
    centered_data[,i] = vec[,i] - mu.a
  }
  for (i in (k.star:N1)){
    centered_data[,i] = vec[,i] - mu.b
  }

  # estimate eigenvalues of covariance operator
  D_mat <- LongRunCovMatrix(centered_data, h = h)
  eigen_struct <- eigen(D_mat, symmetric = TRUE)
  lambda <- Re(eigen_struct$values)/M

  asymp <- function(N) {
    BridgeLam <- matrix(0, M, N)
    for (j in (1:M)) {
      BridgeLam[j, ] <- lambda[j] * (BBridge(0, 0, 0, 1, N - 1)^2)
    }
    max(colSums(BridgeLam))
  }

  Values <- sapply(1:M, function(k) asymp(N1))
  z <- Tn <= Values
  p <- length(z[z == TRUE]) / length(z)

  # Plot
  delta <- mean.a - mean.b
  if (showplot == TRUE) {
    graphics::par(mfrow = c(1, 3))
    graphics::matplot(f, type = "l", col = "grey", main = "Functional Data", ylab = "values")
    for (i in 1:ncol(dat.a)) {
      graphics::lines(dat.a[, i], col = "pink")
    }
    for (i in 1:ncol(dat.b)) {
      graphics::lines(dat.b[, i], col = "lightblue")
    }
    graphics::lines(mean.b, col = "blue")
    graphics::lines(mean.a, col = "red")
    graphics::legend("topleft", c("before", "after"), col = c("blue", "red"), lty = c(1, 1), cex = 0.5)
    graphics::plot(delta, type = "l", main = "Estimated Change Function", ylab = "values")

    graphics::matplot(out$warping_functions, type = "l", col = "grey", main = "Warping Functions", ylab = "values")
    for (i in 1:ncol(warp.b)) {
      graphics::lines(warp.b[, i], col = "pink")
    }
    for (i in 1:ncol(warp.a)) {
      graphics::lines(warp.a[, i], col = "lightblue")
    }
    graphics::lines(warp_mean_a, col = "blue")
    graphics::lines(warp_mean_b, col = "red")
    graphics::legend("topleft", c("before", "after"), col = c("blue", "red"), lty = c(1, 1), cex = 0.5)
  }

  out <- list(
    pvalue = p, change = k.star,
    DataBefore = dat.a, DataAfter = dat.b,
    MeanBefore = mean.a, MeanAfter = mean.b,
    WarpingMeanBefore = warp_mean_a, WarpingMeanAfter = warp_mean_b,
    WarpingBefore = warp.a, WarpingAfter = warp.b,
    change_fun = delta, Sn = Sn, mu=mu
  )

  return(out)
}


#' Elastic Changepoint Detection
#'
#' This function identifies changepoints using a functional PCA
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param pca.method string specifying pca method (options = "combined",
#' "vert", or "horiz", default = "combined")
#' @param pc percentage of cummulation explained variance (default = 0.95)
#' @param d number of monte carlo iterations of Brownian Bridge (default = 1000)
#' @param n_pcs scalar specify number of principal components (default = 5)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param showplot show results plots (default = T)
#' @return Returns a list object containing
#' \item{pvalue}{p value}
#' \item{change}{indice of changepoint}
#' \item{DataBefore}{functions before changepoint}
#' \item{DataAfter}{functions after changepoint}
#' \item{MeanBefore}{mean function before changepoint}
#' \item{MeanAfter}{mean function after changepoint}
#' \item{WarpingBefore}{warping functions before changepoint}
#' \item{WarpingAfter}{warping functions after changepoint}
#' \item{WarpingMeanBefore}{mean warping function before changepoint}
#' \item{WarpingMeanAfter}{mean warping function after changepoint}
#' \item{change_fun}{amplitude change function}
#' \item{Sn}{test statistic values}
#' @keywords srvf alignment changepoint
#' @references J. D. Tucker and D. Yarger, “Elastic Functional Changepoint
#'   Detection of Climate Impacts from Localized Sources”, Envirometrics,
#'   10.1002/env.2826, 2023.
#' @export
elastic_change_fpca <- function(f, time, pca.method = "combined", pc = 0.95, d = 1000, n_pcs = 5, smooth_data=FALSE, sparam=25, showplot = TRUE) {
  M = length(time)
  N1 = dim(f)[2]
  if (smooth_data) f <- smooth.data(f, sparam)

  pca.method1 <- pca.method
  pca.method <- pmatch(pca.method, c("combined", "vert", "horiz")) # 1 - combined, 2 - vert, 3 - horiz
  if (is.na(pca.method)) {
    stop("invalid method selection")
  }

  out <- time_warping(f, time, parallel = T)

  # Calculate PCA -----------------------------------------------------------
  if (pca.method == 1)
    out.pca <- jointFPCA(out, M, showplot=F, C=NULL)
  if (pca.method == 2)
    out.pca <- vertFPCA(out, M, showplot=F)
  if (pca.method == 3)
    out.pca <- horizFPCA(out, M, showplot=F)

  cumm_coef = cumsum(out.pca$latent)/sum(out.pca$latent)
  no = which(cumm_coef >= pc)
  no = no[1]

  lam = 1/out.pca$latent[1:no]
  Sigma = diag(lam)
  eta = out.pca$coef[, 1:no]
  eta_bar = apply(eta,2,sum)

  # compute test statistic
  Sn <- rep(0, N1)
  for (j in (2:N1)) {
    tmp_eta = eta[1:j,]
    tmp = apply(tmp_eta,2,sum)-dim(tmp_eta)[1]*eta_bar

    Sn[j] = 1/N1 * t(tmp) %*% Sigma %*% tmp
  }
  k.star <- min(which(Sn == max(Sn)))
  Tn <- mean(Sn)

  # compute distribution
  Values = rep(0, d)
  for (i in (1:d)){
    B_tmp = matrix(0,no, N1)
    for (j in (1:no)){
      B_tmp[j,] = BBridge(0, 0, 0, 1, N1-1)^2
    }
    Values[i] = mean(apply(B_tmp,2,sum))
  }

  z <- Tn <= Values
  p <- length(z[z == TRUE]) / length(z)
  dat.a <- f[, 1:k.star]
  dat.b <- f[, (k.star + 1):N1]
  warp.a <- out$warping_functions[,1:k.star]
  warp.b <- out$warping_functions[,(k.star + 1):N1]
  mean.a <- rowMeans(out$fn[,1:k.star])
  mean.b <- rowMeans(out$fn[,(k.star + 1):N1])
  out.gam.mean <- SqrtMean(warp.a)
  warp_mean_a = out.gam.mean$gam_mu
  out.gam.mean <- SqrtMean(warp.b)
  warp_mean_b = out.gam.mean$gam_mu
  delta <- mean.a - mean.b

  # Plot
  delta <- mean.a - mean.b
  if (showplot == TRUE) {
    graphics::par(mfrow = c(1, 3))
    graphics::matplot(f, type = "l", col = "grey", main = "Functional Data", ylab = "values")
    for (i in 1:ncol(dat.b)) {
      graphics::lines(dat.b[, i], col = "pink")
    }
    for (i in 1:ncol(dat.a)) {
      graphics::lines(dat.a[, i], col = "lightblue")
    }
    graphics::lines(mean.a, col = "blue")
    graphics::lines(mean.b, col = "red")
    graphics::legend("topleft", c("before", "after"), col = c("blue", "red"), lty = c(1, 1), cex = 0.5)

    graphics::plot(delta, type = "l", main = "Estimated Change Function", ylab = "values")

    graphics::matplot(out$warping_functions, type = "l", col = "grey", main = "Warping Functions", ylab = "values")
    for (i in 1:ncol(warp.b)) {
      graphics::lines(warp.b[, i], col = "pink")
    }
    for (i in 1:ncol(warp.a)) {
      graphics::lines(warp.a[, i], col = "lightblue")
    }
    graphics::lines(warp_mean_a, col = "blue")
    graphics::lines(warp_mean_b, col = "red")
    graphics::legend("topleft", c("before", "after"), col = c("blue", "red"), lty = c(1, 1), cex = 0.5)
  }

  out <- list(
    pvalue = p, change = k.star,
    DataBefore = dat.a, DataAfter = dat.b,
    MeanBefore = mean.a, MeanAfter = mean.b,
    WarpingMeanBefore = warp_mean_a, WarpingMeanAfter = warp_mean_b,
    WarpingBefore = warp.a, WarpingAfter = warp.b,
    change_fun = delta, Sn = Sn
  )

  return(out)
}
