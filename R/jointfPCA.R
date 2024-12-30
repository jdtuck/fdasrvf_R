#' Joint Vertical and Horizontal Functional Principal Component Analysis
#'
#' This function calculates amplitude and phase joint  functional principal component
#' analysis on aligned data
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param no number of principal components to extract (default = 3)
#' @param var_exp compute no based on value percent variance explained (example: 0.95)
#'                will override `no`
#' @param id integration point for f0 (default = midpoint)
#' @param C balance value (default = NULL)
#' @param ci geodesic standard deviations (default = c(-1,0,1))
#' @param showplot show plots of principal directions (default = T)
#' @return Returns a list containing \item{q_pca}{srvf principal directions}
#' \item{f_pca}{f principal directions}
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' \item{mu_psi}{mean psi function}
#' \item{mu_g}{mean g function}
#' \item{id}{point use for f(0)}
#' \item{C}{optimized phase amplitude ratio}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Jung, S. L. a. S. (2016). "Combined Analysis of Amplitude and Phase Variations in Functional Data."
#'        	arXiv:1603.01775.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' jfpca <- jointFPCA(simu_warp, no = 3)
jointFPCA <- function(warp_data,
                      no = 3,
                      var_exp = NULL,
                      id = round(length(warp_data$time) / 2),
                      C = NULL,
                      ci = c(-1, 0, 1),
                      showplot = T) {
  fn <- warp_data$fn
  time <- warp_data$time
  qn <- warp_data$qn
  M <- nrow(qn)
  N <- ncol(qn)
  q0 <- warp_data$q0
  gam <- warp_data$warping_functions
  # Set up for fPCA in q-space
  mq_new <- rowMeans(qn)
  m_new <- sign(fn[id, ]) * sqrt(abs(fn[id, ]))  # scaled version
  mqn <- c(mq_new, mean(m_new))
  qn1 <- rbind(qn, m_new)

  # Calculate Vector Space of warping
  Tgam <- SqrtMean(gam)
  mu_psi <- Tgam$mu
  vec <- Tgam$vec

  # Joint fPCA --------------------------------------------------------------
  jointfPCAd <- function(qn, vec, C = 1, m = 3) {
    M <- nrow(qn)
    N <- ncol(qn)
    time1 <- seq(0, 2, length.out = M + nrow(vec))
    g <- rbind(qn, C * vec)
    mu_q <- rowMeans(qn)

    mu_g <- rowMeans(g)

    K <- stats::cov(t(g))
    out.K <- svd(K, nu = m, nv = m)
    s <- out.K$d
    U <- out.K$u

    a <- matrix(0, N, m)
    for (i in 1:N) {
      for (j in 1:m) {
        a[i, j] <- (g[, i] - mu_g) %*% U[, j]
      }
    }

    qhat <- matrix(mu_q, M, N) + U[1:M, 1:m] %*% t(a)
    vechat <- U[(M + 1):length(time1), 1:m] %*% t(a / C)
    psihat <- apply(vechat, 2, function(x, mu_psi) {
      exp_map(mu_psi, x)
    }, mu_psi)
    gamhat <- apply(psihat, 2, function(x, time) {
      gam_mu = cumtrapz(time, x * x)
      gam_mu = (gam_mu - min(gam_mu)) / (max(gam_mu) - min(gam_mu))
    }, seq(0, 1, length.out = M - 1))

    return(list(
      qhat = qhat,
      gamhat = gamhat,
      a = a,
      U = U[, 1:m],
      s = s[1:m],
      eigs = s,
      mu_g = mu_g,
      cov = K,
      g = g
    ))
  }


  # Find C ------------------------------------------------------------------
  findC <- function(C, qn, vec, q0, m) {
    out.pca <- jointfPCAd(qn, vec, C, m)
    M <- nrow(qn)
    N <- ncol(qn)
    time <- seq(0, 1, length.out = M - 1)

    d <- rep(0, N)
    for (i in 1:N) {
      tmp <- warp_q_gamma(out.pca$qhat[1:(M - 1), i], time, invertGamma(out.pca$gamhat[, i]))
      d[i] <- sum(trapz(time, (tmp - q0[, i]) ^ 2))
    }

    return(sum(d ^ 2) / N)
  }

  m <- M
  if (is.null(C))
    C <- stats::optimize(
      findC,
      c(0, 1e4),
      qn = qn1,
      vec = vec,
      q0 = q0,
      m = m
    )$minimum

  # Final PCA ---------------------------------------------------------------
  out.pca <- jointfPCAd(qn1, vec, C, m = m)

  if (!is.null(var_exp)) {
    cumm_coef <- cumsum(out.pca$s) / sum(out.pca$s)
    tmp = which(cumm_coef <= var_exp)
    no = tmp[length(tmp)]
  }

  # geodesic paths
  q_pca <- array(0, dim = c(M, length(ci), no))
  f_pca <- array(0, dim = c(M, length(ci), no))
  N1 <- nrow(out.pca$U)
  for (j in 1:no) {
    for (i in 1:length(ci)) {
      qhat <- mqn + out.pca$U[1:(M + 1), j] * ci[i] * sqrt(out.pca$s[j])
      vechat <- out.pca$U[(M + 2):N1, j] * (ci[i] * sqrt(out.pca$s[j])) /
        C
      psihat <- exp_map(mu_psi, vechat)
      gamhat <- cumtrapz(seq(0, 1, length.out = M), psihat * psihat)
      gamhat = (gamhat - min(gamhat)) / (max(gamhat) - min(gamhat))
      if (sum(vechat) == 0)
        gamhat <- seq(0, 1, length.out = M)
      if (id == 1)
        fhat <- cumtrapz(time, qhat[1:M] * abs(qhat[1:M])) + sign(qhat[M +
                                                                         1]) * (qhat[M + 1] ^ 2)
      else
        fhat <- cumtrapzmid(time, qhat[1:M] * abs(qhat[1:M]), sign(qhat[M +
                                                                          1]) * (qhat[M + 1] ^ 2), id)
      f_pca[, i, j] <- warp_f_gamma(fhat, seq(0, 1, length.out = M), gamhat)
      q_pca[, i, j] <- warp_q_gamma(qhat[1:M], seq(0, 1, length.out =
                                                     M), gamhat)
    }
  }

  jfpca <- list()
  jfpca$q_pca <- q_pca
  jfpca$f_pca <- f_pca
  jfpca$latent <- out.pca$s[1:no]
  jfpca$eigs <-  out.pca$eigs
  jfpca$coef <- out.pca$a[, 1:no]
  jfpca$U <- out.pca$U[, 1:no]
  jfpca$mu_psi <- mu_psi
  jfpca$mu_g <- out.pca$mu_g
  jfpca$id <- id
  jfpca$C <- C
  jfpca$stds <- ci
  jfpca$time <- time
  jfpca$g <- out.pca$g
  jfpca$cov <- out.pca$cov
  jfpca$warp_data <- warp_data

  class(jfpca) <- "jfpca"

  if (showplot) {
    plot(jfpca)
  }

  return(jfpca)
}


#' Joint Vertical and Horizontal Functional Principal Component Analysis
#'
#' This function calculates amplitude and phase joint functional principal component
#' analysis on aligned data using the SRVF framework using MFPCA and h representation
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param var_exp compute no based on value percent variance explained (default: 0.99)
#'                will override `no`
#' @param id integration point for f0 (default = midpoint)
#' @param C balance value (default = NULL)
#' @param ci geodesic standard deviations (default = c(-1,0,1))
#' @param srvf use srvf (default = TRUE)
#' @param showplot show plots of principal directions (default = T)
#' @return Returns a list containing \item{q_pca}{srvf principal directions}
#' \item{f_pca}{f principal directions}
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' \item{mu_psi}{mean psi function}
#' \item{mu_g}{mean g function}
#' \item{id}{point use for f(0)}
#' \item{C}{optimized phase amplitude ratio}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Jung, S. L. a. S. (2016). "Combined Analysis of Amplitude and Phase Variations in Functional Data."
#'        	arXiv:1603.01775.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' jfpcah <- jointFPCAh(simu_warp)
jointFPCAh <- function(warp_data,
                       var_exp = 0.99,
                       id = round(length(warp_data$time) / 2),
                       C = NULL,
                       ci = c(-1, 0, 1),
                       srvf = TRUE,
                       showplot = TRUE) {
  fn <- warp_data$fn
  time <- warp_data$time
  qn <- warp_data$qn
  M <- nrow(qn)
  N <- ncol(qn)
  q0 <- warp_data$q0
  gam <- warp_data$warping_functions

  # Set up for fPCA in q-space
  if (srvf) {
    mq_new <- rowMeans(qn)
    m_new <- sign(fn[id, ]) * sqrt(abs(fn[id, ]))  # scaled version
    mqn <- c(mq_new, mean(m_new))
    qn1 <- rbind(qn, m_new)
  } else {
    mqn <- rowMeans(fn)
    q0 <- warp_data$f0
    qn1 <- fn
  }


  # Calculate Vector Space of warping
  h <- gam_to_h(gam, smooth=FALSE)

  # Joint fPCA --------------------------------------------------------------
  jointfPCAhd <- function(qn, h, C = 1, var_exp = NULL) {
    M <- nrow(qn)
    N <- ncol(qn)

    # Run Univariate fPCA
    # q space
    mqn <- rowMeans(qn)

    K <- stats::cov(t(qn))
    out.K <- svd(K)
    s <- out.K$d
    U <- out.K$u

    cumm_coef <- cumsum(s) / sum(s)
    tmp = which(cumm_coef <= var_exp)
    no_q = tmp[length(tmp)]

    c.o = t(qn - mqn) %*% U
    c = c.o[, 1:no_q, drop = F]
    U = U[, 1:no_q, drop = F]

    if (no_q == 1){
      U = matrix(U)
      c = matrix(c)
    }

    # h space
    hc = C * h
    mh = rowMeans(hc)
    Kh = stats::cov(t(hc))

    out.Kh <- svd(Kh)
    sh <- out.Kh$d
    Uh <- out.Kh$u

    cumm_coef <- cumsum(sh) / sum(sh)
    tmp = which(cumm_coef <= var_exp)
    no_h = tmp[length(tmp)]

    ch.o = t(hc - mh) %*% Uh
    ch = ch.o[, 1:no_h, drop = F]
    Uh = Uh[, 1:no_h, drop = F]

    # Run Multivariate fPCA
    Xi = cbind(c, ch)
    Xi.o = cbind(c.o, ch.o)
    Z = 1 / (nrow(Xi) - 1) * t(Xi) %*% Xi

    out.Z <- svd(Z)
    sz <- out.Z$d
    Uz <- out.Z$u

    Z.o = 1 / (nrow(Xi.o) - 1) * t(Xi.o) %*% Xi.o
    out.Zo <- svd(Z.o)

    cz = Xi %*% Uz

    if (no_q == 1){
      Psi_q = U %*% t(matrix(Uz[1:no_q, ]))
    } else {
      Psi_q = U %*% Uz[1:no_q, ]
    }

    Psi_h = Uh %*% Uz[(no_q + 1):nrow(Uz), ]

    hhat = Psi_h %*% t(cz)
    gamhat = h_to_gam(hhat / C)
    qhat = Psi_q %*% t(cz) + mqn

    return(
      list(
        qhat = qhat,
        gamhat = gamhat,
        cz = cz,
        Psi_q = Psi_q,
        Psi_h = Psi_h,
        sz = sz,
        eigs = out.Zo$d,
        U = U,
        Uh = Uh,
        Uz = Uz
      )
    )
  }


  # Find C ------------------------------------------------------------------
  findCh <- function(C, qn, h, q0, var_exp, srvf) {
    out.pca <- jointfPCAhd(qn, h, C, var_exp)
    M <- nrow(qn)
    N <- ncol(qn)

    d <- rep(0, N)
    for (i in 1:N) {
      if (srvf) {
        time <- seq(0, 1, length.out = M - 1)
        tmp <- warp_q_gamma(out.pca$qhat[1:(M - 1), i], time, invertGamma(out.pca$gamhat[, i]))
      } else {
        time <- seq(0, 1, length.out = M)
        tmp <- warp_f_gamma(out.pca$qhat[, i], time, invertGamma(out.pca$gamhat[, i]))
      }

      d[i] <- sum(trapz(time, (tmp - q0[, i]) ^ 2))
    }

    return(mean(d))
  }

  findCh(.1,qn1, h, q0, 0.99, FALSE)
  if (is.null(C))
    C <- stats::optimize(
      findCh,
      c(0, 1e4),
      qn = qn1,
      h = h,
      q0 = q0,
      var_exp = 0.99,
      srvf = srvf
    )$minimum

  # Final PCA ---------------------------------------------------------------
  out.pca <- jointfPCAhd(qn1, h, C, var_exp = var_exp)

  hc <- C * h
  mh = rowMeans(hc)

  # geodesic paths
  no = ncol(out.pca$cz)
  q_pca <- array(0, dim = c(M, length(ci), no))
  f_pca <- array(0, dim = c(M, length(ci), no))
  N1 <- nrow(out.pca$U)
  for (j in 1:no) {
    for (i in 1:length(ci)) {
      qhat <- mqn + out.pca$Psi_q[, j] * (ci[i] * sqrt(out.pca$sz[j]))
      hhat = out.pca$Psi_h[, j] * (ci[i] * sqrt(out.pca$sz[j])) / C
      gamhat = h_to_gam(hhat)
      gamhat = (gamhat - min(gamhat)) / (max(gamhat) - min(gamhat))

      if (srvf) {
        if (id == 1)
          fhat <- cumtrapz(time, qhat[1:M] * abs(qhat[1:M])) + sign(qhat[M +
                                                                           1]) * (qhat[M + 1] ^ 2)
        else
          fhat <- cumtrapzmid(time, qhat[1:M] * abs(qhat[1:M]), sign(qhat[M +
                                                                            1]) * (qhat[M + 1] ^ 2), id)
        f_pca[, i, j] <- warp_f_gamma(fhat, seq(0, 1, length.out = M), gamhat)
        q_pca[, i, j] <- warp_q_gamma(qhat[1:M], seq(0, 1, length.out =
                                                       M), gamhat)
      } else {
        f_pca[, i, j] <- warp_f_gamma(qhat, seq(0, 1, length.out = M), gamhat)
        q_pca[, i, j] <- f_to_srvf(f_pca[, i, j], seq(0, 1, length.out = M))
      }

    }
  }

  jfpcah <- list()
  jfpcah$q_pca <- q_pca
  jfpcah$f_pca <- f_pca
  jfpcah$eigs <- out.pca$eigs
  jfpcah$latent <- out.pca$sz[1:no]
  jfpcah$coef <- out.pca$cz[, 1:no]
  jfpcah$U_q <- out.pca$Psi_q
  jfpcah$U_h <- out.pca$Psi_h
  jfpcah$id <- id
  jfpcah$C <- C
  jfpcah$stds <- ci
  jfpcah$time <- time
  jfpcah$h <- h
  jfpcah$qn1 <- qn1
  jfpcah$mqn <- mqn
  jfpcah$U <- out.pca$U
  jfpcah$U1 <- out.pca$Uh
  jfpcah$Uz <- out.pca$Uz
  jfpcah$mh <- mh
  jfpcah$warp_data <- warp_data

  class(jfpcah) <- "jfpcah"

  if (showplot) {
    plot(jfpcah)
  }

  return(jfpcah)
}
