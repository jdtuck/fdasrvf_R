#' Group-wise function alignment and PCA Extractions
#'
#' This function aligns a collection of functions while extracting principal
#' components.
#'
#' @param f A numeric matrix of shape \eqn{M \times N} specifying a sample of
#'   \eqn{N} \eqn{1}-dimensional curves observed on a grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   functions `f` have been evaluated.
#' @param num_comp An integer value specifying the number of principal
#'   components to extract. Defaults to `3L`.
#' @param showplot A boolean specifying whether to display plots along the way.
#'   Defaults to `TRUE`.
#' @param smooth_data A boolean specifying whether to smooth data using box
#'   filter. Defaults to `FALSE`.
#' @param sparam An integer value specifying the number of times to apply box
#'   filter. Defaults to `25L`. This argument is only used if `smooth_data ==
#'   TRUE`.
#' @param parallel A boolean specifying whether computations should run in
#'   parallel. Defaults to `FALSE`.
#' @param cores An integer value specifying the number of cores to use for
#'   parallel computations. Defaults to `NULL` in which case it uses all
#'   available cores but one. This argument is only used when `parallel ==
#'   TRUE`.
#' @param max_iter An integer value specifying the maximum number of iterations.
#'   Defaults to `51L`.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#'
#' @return A list with the following components:
#'
#' - `f0`: A numeric matrix of shape \eqn{M \times N} storing the original
#' functions;
#' - `fn`: A numeric matrix of the same shape as `f0` storing the aligned
#' functions;
#' - `qn`: A numeric matrix of the same shape as `f0` storing the aligned
#' SRSFs;
#' - `q0`: A numeric matrix of the same shape as `f0` storing the SRSFs of the
#' original functions;
#' - `mqn`: A numeric vector of length \eqn{M} storing the mean SRSF;
#' - `gam`: A numeric matrix of the same shape as `f0` storing the estimated
#' warping functions;
#' - `vfpca`: A list storing information about the vertical PCA with the
#' following components:
#'
#'   - `q_pca`: A numeric matrix of shape \eqn{(M + 1) \times 5 \times
#'   \mathrm{num\_comp}} storing the first \eqn{3} principal directions in SRSF
#'   space; the first dimension is \eqn{M + 1} because, in SRSF space, the
#'   original functions are represented by the SRSF and the initial value of the
#'   functions.
#'   - `f_pca`: A numeric matrix of shape \eqn{M \times 5 \times
#'   \mathrm{num\_comp}} storing the first \eqn{3} principal directions in
#'   original space;
#'   - `latent`: A numeric vector of length \eqn{M + 1} storing the singular
#'   values of the SVD decomposition in SRSF space;
#'   - `coef`: A numeric matrix of shape \eqn{N \times \mathrm{num\_comp}}
#'   storing the scores of the \eqn{N} original functions on the first
#'   `num_comp` principal components;
#'   - `U`: A numeric matrix of shape \eqn{(M + 1) \times (M + 1)} storing the
#'   eigenvectors associated with the SVD decomposition in SRSF space.
#'
#' - `Dx`: A numeric vector of length `max_iter` storing the value of the cost
#' function at each iteration.
#'
#' @keywords pca
#' @concept srvf alignment
#'
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using phase and amplitude separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#'
#' @export
#' @examples
#' \dontrun{
#'   out <- align_fPCA(simu_data$f, simu_data$time)
#' }
align_fPCA <- function(f, time,
                       num_comp = 3L,
                       showplot = TRUE,
                       smooth_data = FALSE,
                       sparam = 25L,
                       parallel = FALSE,
                       cores = NULL,
                       max_iter = 51L,
                       lambda = 0.0) {
  if (parallel) {
    if (is.null(cores))
      cores <- max(parallel::detectCores() - 1, 1)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  } else
    foreach::registerDoSEQ()

  cli::cli_alert_info("Initializing...")

  binsize <- mean(diff(time))
  eps <- .Machine$double.eps
  M <- nrow(f)
  N <- ncol(f)
  f0 <- f
  coef <- -2:2
  NP <- 1:num_comp  # number of principal components
  Nstd <- length(coef)

  if (smooth_data) f <- smooth.data(f, sparam = sparam)

  if (showplot) {
    graphics::matplot(time, f, type = "l")
    graphics::title(main = "Original data")
  }

  # Compute q-function of the plot
  tmp <- gradient.spline(f, binsize)
  f <- tmp$f
  q <- tmp$g / sqrt(abs(tmp$g) + eps)

  # PCA Step
  mnq <- rowMeans(q)
  dqq <- sqrt(colSums((q - matrix(mnq, nrow = M, ncol = N))^2))
  min_ind <- which.min(dqq)
  mq <- q[, min_ind]
  qhat_cent <- q - matrix(mq, nrow = M, ncol = N)
  K <- 1/M * qhat_cent %*% t(qhat_cent)
  out <- svd2(K)
  s <- out$d
  U <- out$u

  alpha_i <- matrix(0, nrow = num_comp, ncol = N)
  for (ii in 1:num_comp) {
    for (jj in 1:N)
      alpha_i[ii, jj] <- simpson(time, qhat_cent[, jj] * U[, ii])
  }

  tmp <- U[, 1:num_comp] %*% alpha_i
  if (smooth_data) {
    mq <- mq / pvecnorm(mq, 2)
    for (ii in 1:N) {
      if (sum(tmp[, ii]) != 0)
        tmp[, ii] <- tmp[, ii] / pvecnorm(tmp[, ii], 2)
    }
  }
  qhat <- matrix(mq, nrow = M, ncol = N) + tmp

  # Matching Step
  gam_k <- foreach::foreach(k = 1:N, .combine = cbind, .packages="fdasrvf") %dopar% {
    gam0 <- optimum.reparam(qhat[, k], time, q[, k], time, lambda = lambda)
  }

  cli::cli_alert_info("Aligning {N} function{?s} in SRVF space to {num_comp} fPCA component{?s}...")

  tmp <- matrix(0, nrow = M, ncol = max_iter + 2)
  tmp[, 1] <- mq
  mq <- tmp
  tmp <- array(0, dim = c(M, N, max_iter + 2))
  tmp[, , 1] <- f
  f <- tmp
  tmp <- array(0, dim = c(M, N, max_iter + 2))
  tmp[, , 1] <- q
  q <- tmp
  tmp <- array(0, dim = c(M, N, max_iter + 2))
  tmp[, , 1] <- gam_k
  gam <- tmp

  Dx <- rep(0, max_iter)

  psi <- sqrt(diff(gam_k) * (M - 1))
  Dx1 <- rep(0, N)
  for (ii in 1:N) {
    acos_input <- mean(psi[, ii])
    if (acos_input >  1) acos_input <-  1
    if (acos_input < -1) acos_input <- -1
    Dx1[ii] <- Re(acos(acos_input))
  }
  Dx[1] <- max(Dx1)

  for (r in 2:max_iter) {
    cli::cli_alert_info("Starting iteration {r - 1}...")

    if (r == max_iter)
      cli::cli_alert_info("The maximal number of iterations has been reached.")

    # PCA Step
    f_temp <- matrix(0, nrow = M, ncol = N)
    q_temp <- matrix(0, nrow = M, ncol = N)

    for (k in 1:N) {
      gam_k <- gam[, k, r - 1]
      f_temp[,k] <- stats::approx(
        x = time,
        y = f[, k, r - 1],
        xout = (time[M] - time[1]) * gam_k + time[1]
      )$y
      q_temp[, k] <- f_to_srvf(f_temp[, k], time)
    }

    f[, , r] <- f_temp
    q[, , r] <- q_temp
    mq[, r] <- rowMeans(q_temp)

    K <- stats::cov(t(q[, , r]))
    out <- svd(K)
    s <- out$d
    U <- out$u

    qhat_cent <- scale(t(q[, , r]), scale = FALSE)
    qhat_cent <- t(qhat_cent)
    alpha_i <- matrix(0, nrow = num_comp, ncol = N)

    for (ii in 1:num_comp) {
      for (jj in 1:N)
        alpha_i[ii, jj] <- simpson(time, qhat_cent[, jj] * U[, ii])
    }

    tmp <- U[, 1:num_comp] %*% alpha_i
    mq_c <- mq[, r]
    if (smooth_data) {
      mq_c <- mq[, r] / pvecnorm(mq[, r], 2)
      for (ii in 1:N) {
        if (sum(tmp[, ii]) != 0)
          tmp[, ii] <- tmp[, ii] / pvecnorm(tmp[, ii], 2)
      }
    }
    qhat <- matrix(mq_c, nrow = M, ncol = N) + tmp

    # Matching Step
    gam_k <- foreach::foreach(k = 1:N, .combine = cbind, .packages = "fdasrvf") %dopar% {
      gam0 <- optimum.reparam(qhat[, k], time, q[, k, r], time, lambda = lambda)
    }
    gam[, , r] <- gam_k

    psi <- sqrt(diff(gam_k) * (M - 1))
    Dx1 <- rep(0, N)
    for (ii in 1:N) {
      acos_input <- mean(psi[, ii])
      if (acos_input >  1) acos_input <-  1
      if (acos_input < -1) acos_input <- -1
      Dx1[ii] <- Re(acos(acos_input))
    }
    Dx[r] <- max(Dx1)

    if (abs(Dx[r] - Dx[r - 1]) < 1e-4 || r >= max_iter)
      break
  }

  # Aligned data & stats
  fn <- f[, , r]
  qn <- q[, , r]
  q0 <- q[, , 1]
  mean_f0 <- rowMeans(f[, , 1])
  std_f0 = apply(f[, , 1], 1, stats::sd)
  mqn <- mq[, r]
  gamf <- gam[, , 1]
  for (k in 2:r) {
    gam_k <- gam[, , k]
    for (l in 1:N)
      gamf[, l] <- stats::approx(
        x = time,
        y = gamf[, l],
        xout = (time[M] - time[1]) * gam_k[, l] + time[1]
      )$y
  }

  # Center Mean
  gamI <- SqrtMeanInverse(gamf)
  gamI_dev <- gradient(gamI, 1 / (M - 1))
  mqn <- stats::approx(
    x = time,
    y = mqn,
    xout = (time[M] - time[1]) * gamI + time[1]
  )$y * sqrt(gamI_dev)

  for (k in 1:N) {
    qn[, k] <- stats::approx(
      x = time,
      y = qn[, k],
      xout = (time[M] - time[1]) * gamI + time[1]
    )$y * sqrt(gamI_dev)
    fn[, k] <- stats::approx(
      x = time,
      y = fn[, k],
      xout = (time[M] - time[1]) * gamI + time[1]
    )$y
    gamf[, k] <- stats::approx(
      x = time,
      y = gamf[, k],
      xout = (time[M] - time[1]) * gamI + time[1]
    )$y
  }

  mean_fn <- rowMeans(fn)
  std_fn <- apply(fn, 1, stats::sd)

  # Get Final PCA
  mq_new <- rowMeans(qn)
  id <- round(M / 2)
  m_new <- sign(fn[id, ]) * sqrt(abs(fn[id, ])) # scaled version
  mqn2 <- c(mq_new, mean(m_new))
  K <- stats::cov(t(rbind(qn, m_new)))
  out <- svd(K)
  s <- out$d
  stdS <- sqrt(s)
  U <- out$u

  # compute the PCA in the q domain
  q_pca <- array(0, dim = c(M + 1, Nstd, num_comp))
  for (k in NP) {
    for (i in 1:Nstd)
      q_pca[, i, k] <- mqn2 + coef[i] * stdS[k] * U[, k]
  }

  # compute the correspondence to the original function domain
  f_pca <- array(0, dim = c(M, Nstd, num_comp))
  for (k in NP) {
    for (i in 1:Nstd)
      f_pca[, i, k] <- cumtrapzmid(
        x = time,
        y = q_pca[1:M, i, k] * abs(q_pca[1:M, i, k]),
        c = sign(q_pca[M + 1, i, k]) * q_pca[M + 1, i, k]^2,
        mid = id
      )
  }

  c <- matrix(0, N, num_comp)
  for (k in NP) {
    for (i in 1:N)
      c[i, k] <- sum((c(qn[, i], m_new[i]) - mqn2) * U[, k])
  }

  vfpca = list()
  vfpca$q_pca <- q_pca
  vfpca$f_pca <- f_pca
  vfpca$latent <- s
  vfpca$coef <- c
  vfpca$U <- U

  if (showplot) {
    graphics::matplot(
      x = (0:(M - 1)) / (M - 1),
      y = gamf,
      type = "l",
      main = "Warping functions",
      xlab = "Time"
    )

    graphics::matplot(
      x = time,
      y = fn,
      type = "l",
      main = bquote(paste("Warped Data ",lambda == .(lambda)))
    )

    graphics::matplot(
      x = time,
      y = cbind(mean_f0, mean_f0 + std_f0, mean_f0 - std_f0),
      type = "l",
      lty = 1,
      col = c("blue", "red", "green"),
      ylab = "",
      main = bquote(paste("Original Data: ", Mean %+-% STD))
    )
    graphics::legend(
      x = "topright",
      inset = 0.01,
      legend = c("Mean", "Mean + STD", "Mean - STD"),
      col = c("blue", "red", "green"),
      lty = 1
    )

    graphics::matplot(
      x = time,
      y = cbind(mean_fn, mean_fn + std_fn, mean_fn - std_fn),
      type = "l",
      lty = 1,
      col = c("blue", "red", "green"),
      ylab = "",
      main = bquote(paste(
        "Warped Data: ",
        lambda == .(lambda),
        ": ",
        Mean %+-% STD
      ))
    )
    graphics::legend(
      x = "topright",
      inset = 0.01,
      legend = c("Mean", "Mean + STD", "Mean - STD"),
      col = c("blue", "red", "green"),
      lty = 1
    )

    if (num_comp > 3) num_comp <- 3
    graphics::layout(matrix(
      data = 1:(2 * num_comp),
      nrow = 2,
      ncol = num_comp,
      byrow = TRUE
    ))
    dims <- dim(q_pca)
    for (ii in 1:num_comp) {
      graphics::matplot(time, vfpca$q_pca[1:M, , ii], type = "l")
      graphics::title(main = sprintf("q domain: PD %d", ii))
    }
    for (ii in 1:num_comp) {
      graphics::matplot(time, vfpca$f_pca[, , ii], type = "l")
      graphics::title(main = sprintf("f domain: PD %d", ii))
    }

    graphics::layout(1)
    cumm_coef <- 100 * cumsum(vfpca$latent) / sum(vfpca$latent)
    plot(
      x = cumm_coef,
      type = "l",
      col = "blue",
      main = "Coefficient Cumulative Percentage",
      ylab = "Percentage"
    )
  }

  if (parallel) parallel::stopCluster(cl)

  list(
    f0 = f[, , 1],
    q0 = q0,
    fn = fn,
    qn = qn,
    mqn = mqn,
    gam = gamf,
    vfpca = vfpca,
    Dx = Dx
  )
}
