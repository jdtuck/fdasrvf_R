#' Alignment of univariate functional data
#'
#' This function aligns a collection of \eqn{1}-dimensional curves that are
#' observed on the same grid.
#'
#' @param f A numeric matrix of shape \eqn{M \times N} specifying a sample of
#'   \eqn{N} curves observed on a grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the common grid on
#'   which all curves `f` have been observed.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param penalty_method A string specifying the penalty term used in the
#'   formulation of the cost function to minimize for alignment. Choices are
#'   `"roughness"` which uses the norm of the second derivative, `"geodesic"`
#'   which uses the geodesic distance to the identity and `"norm"` which uses
#'   the Euclidean distance to the identity. Defaults to `"roughness"`.
#' @param centroid_type A string specifying the type of centroid to align to.
#'   Choices are `"mean"` or `"median"`. Defaults to `"mean"`.
#' @param center_warpings A boolean specifying whether to center the estimated
#'   warping functions. Defaults to `TRUE`.
#' @param smooth_data A boolean specifying whether to smooth curves using a box
#'   filter. Defaults to `FALSE`.
#' @param sparam An integer value specifying the number of times to apply the
#'   box filter. Defaults to `25L`. This is used only when `smooth_data = TRUE`.
#' @param parallel A boolean specifying whether to run calculations in parallel.
#'   Defaults to `FALSE`.
#' @param optim_method A string specifying the algorithm used for optimization.
#'   Choices are `"DP"`, `"DP2"` and `"RBFGS"`. Defaults to `"DP"`.
#' @param max_iter An integer value specifying the maximum number of iterations.
#'   Defaults to `20L`.
#'
#' @return An object of class `fdawarp` which is a list with the following
#'   components:
#'
#' - `time`: a numeric vector of length \eqn{M} storing the original grid;
#' - `f0`: a numeric matrix of shape \eqn{M \times N} storing the original
#' sample of \eqn{N} functions observed on a grid of size \eqn{M};
#' - `q0`: a numeric matrix of the same shape as `f0` storing the original
#' SRSFs;
#' - `fn`: a numeric matrix of the same shape as `f0` storing the aligned
#' functions;
#' - `qn`: a numeric matrix of the same shape as `f0` storing the aligned SRSFs;
#' - `fmean`: a numeric vector of length \eqn{M} storing the mean or median
#' curve;
#' - `mqn`: a numeric vector of length \eqn{M} storing the mean or median SRSF;
#' - `warping_functions`: a numeric matrix of the same shape as `f0` storing the
#' estimated warping functions;
#' - `original_variance`: a numeric value storing the variance of the original
#' sample;
#' - `amplitude_variance`: a numeric value storing the variance in amplitude of
#' the aligned sample;
#' - `phase_variance`: a numeric value storing the variance in phase of the
#' aligned sample;
#' - `qun`: a numeric vector of maximum length `max_iter + 2` storing the values
#' of the cost function after each iteration;
#' - `lambda`: the input parameter `lambda` which specifies the elasticity;
#' - `centroid_type`: the input centroid type;
#' - `optim_method`: the input optimization method;
#' - `inverse_average_warping_function`: the inverse of the mean estimated
#' warping function;
#' - `rsamps`: TO DO.
#'
#' @keywords srsf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using Fisher-Rao metric,
#'   arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using phase and amplitude Separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#'
#' @export
#' @examples
#' \dontrun{
#'   out <- time_warping(simu_data$f, simu_data$time)
#' }
time_warping <- function(f, time,
                         lambda = 0.0,
                         penalty_method = c("roughness", "geodesic", "norm"),
                         centroid_type = c("mean", "median"),
                         center_warpings = TRUE,
                         smooth_data = FALSE,
                         sparam = 25L,
                         parallel = FALSE,
                         optim_method = c("DP", "DP2", "RBFGS"),
                         max_iter = 20L) {
  penalty_method <- rlang::arg_match(penalty_method)
  centroid_type <- rlang::arg_match(centroid_type)
  optim_method <- rlang::arg_match(optim_method)

  if (parallel) {
    cores <- max(parallel::detectCores() - 1, 1)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  } else
    foreach::registerDoSEQ()

  cli::cli_alert_info("Using lambda = {lambda}")

  binsize <- mean(diff(time))
  eps <- .Machine$double.eps
  M <- nrow(f)
  N <- ncol(f)
  f0 <- f
  w <- 0.0

  if (smooth_data) f <- smooth.data(f, sparam)

  # Compute q-function of the functional data
  tmp <- gradient.spline(f, binsize, smooth_data)
  f <- tmp$f
  q <- tmp$g / sqrt(abs(tmp$g) + eps)

  cli::cli_alert_info("Initializing...")

  mnq <- rowMeans(q)
  dqq <- sqrt(colSums((q - matrix(mnq, ncol = N, nrow = M))^2))
  min_ind <- which.min(dqq)
  mq <- q[, min_ind]
  mf <- f[, min_ind]

  gam <- foreach::foreach(n = 1:N, .combine = cbind, .packages = "fdasrvf") %dopar% {
    gam_tmp <- optimum.reparam(
      Q1 = mq, T1 = time, Q2 = q[, n], T2 = time,
      lambda = lambda, pen = penalty_method, method = optim_method, w = w,
      f1o = mf[1], f2o = f[1, n]
    )
  }

  gamI <- SqrtMeanInverse(gam)
  gamI_dev <- gradient(gamI, 1 / (M - 1))
  gam <- t(gam)

  mf <- stats::approx(time, mf, xout = (time[M] - time[1]) * gamI + time[1])$y
  mq <- f_to_srvf(mf, time)
  mq[is.nan(mq)] <- 0

  cli::cli_alert_info("Computing Karcher {centroid_type} of {N} functions in SRSF space...")

  ds <- rep(0, max_iter + 2)
  ds[1] <- Inf

  tmp <- matrix(0, nrow = M, ncol = max_iter + 2)
  tmp[, 1] <- mq
  mq <- tmp

  tmp <- matrix(0, nrow = M, ncol = max_iter + 2)
  tmp[, 1] <- mf
  mf <- tmp

  tmp <- array(0, dim = c(M, N, max_iter + 2))
  tmp[, , 1] <- f
  f <- tmp

  tmp <- array(0, dim = c(M, N, max_iter + 2))
  tmp[, , 1] <- q
  q <- tmp

  qun <- rep(0, max_iter + 2)
  qun[1] <- pvecnorm(mq[, 1] - q[, min_ind, 1], 2) / pvecnorm(q[, min_ind, 1], 2)
  stp <- .3

  for (r in 1:max_iter) {
    cli::cli_alert_info("Entering iteration {r}...")

    if (r == max_iter)
      cli::cli_alert_warning("The maximal number of iterations has been reached.")

    # Matching Step
    outfor <- foreach::foreach(n = 1:N, .combine = cbind, .packages='fdasrvf') %dopar% {
      gam <- optimum.reparam(
        Q1 = mq[, r], T1 = time, Q2 = q[, n, 1], T2 = time,
        lambda = lambda, pen = penalty_method, method = optim_method, w = w,
        f1o = mf[1, r], f2o = f[1, n, 1]
      )
      gam_dev <- gradient(gam, 1 / (M - 1))
      f_temp <- stats::approx(
        time, f[, n, 1],
        xout = (time[M] - time[1]) * gam + time[1]
      )$y
      q_temp <- f_to_srvf(f_temp, time)
      v <- q_temp - mq[, r]
      d <- sqrt(trapz(time, v * v))
      vtil <- v / d
      dtil <- 1 / d

      list(gam, gam_dev, q_temp, f_temp, vtil, dtil)
    }

    gam <- unlist(outfor[1, ])
    dim(gam) <- c(M, N)
    gam <- t(gam)
    gam_dev <- unlist(outfor[2, ])
    dim(gam_dev) <- c(M, N)
    gam_dev <- t(gam_dev)
    q_temp <- unlist(outfor[3, ])
    dim(q_temp) <- c(M, N)
    f_temp <- unlist(outfor[4, ])
    dim(f_temp) <- c(M, N)
    q[, , r + 1] <- q_temp
    f[, , r + 1] <- f_temp
    tmp <- (1 - sqrt(gam_dev)) ^ 2
    vtil <- unlist(outfor[5, ])
    dim(vtil) <- c(M, N)
    dtil <- unlist(outfor[6, ])
    dim(dtil) <- c(1, N)

    if (centroid_type == "mean") {
      ds_tmp <- sum(trapz(time, (matrix(mq[, r], M, N) - q[, , r + 1]) ^ 2)) +
        lambda * sum(trapz(time, t(tmp)))

      if (is.complex(ds_tmp))
        ds[r + 1] <- abs(ds_tmp)
      else
        ds[r + 1] <- ds_tmp

      # Minimization Step
      # compute the mean of the matched function
      mq[, r + 1] <- rowMeans(q[, , r + 1])
      mf[, r + 1] <- rowMeans(f[, , r + 1])
    }

    if (centroid_type == "median") {
      ds_tmp <- sqrt(sum(trapz(time, (q[, , r + 1] - matrix(mq[, r], M, N))^2))) +
        lambda * sum(trapz(time, t(tmp)))

      if (is.complex(ds_tmp))
        ds[r + 1] <- abs(ds_tmp)
      else
        ds[r + 1] <- ds_tmp

      # Minimization Step
      # compute the median of the matched function
      vbar <- rowSums(vtil) * sum(dtil) ^ (-1)
      mq[, r + 1] <- mq[, r] + stp * vbar
      mf[, r + 1] <- stats::median(f[1, , 1]) +
        cumtrapz(time, mq[, r + 1] * abs(mq[, r + 1]))
    }

    qun[r + 1] <- pvecnorm(mq[, r + 1] - mq[, r], 2) / pvecnorm(mq[, r], 2)

    if (qun[r + 1] - qun[r] < 1.0e-4 * qun[r])
      break
  }

  # One last run, centering of gam
  gamI <- SqrtMeanInverse(t(gam))
  gamI_dev <- gradient(gamI, 1 / (M - 1))

  if (center_warpings) {
    r <- r + 1
    outfor <- foreach::foreach(n = 1:N, .combine = cbind, .packages = "fdasrvf") %dopar% {
      gam <- optimum.reparam(
        Q1 = mq[, r], T1 = time, Q2 = q[, n, 1], T2 = time,
        lambda = lambda, pen = penalty_method, method = optim_method, w = w,
        f1o = mf[1, r], f2o = f[1, n, 1]
      )
      gam_dev <- gradient(gam, 1 / (M - 1))
      list(gam, gam_dev)
    }

    gam <- unlist(outfor[1, ])
    dim(gam) <- c(M, N)
    gam <- t(gam)
    gam_dev <- unlist(outfor[2, ])
    dim(gam_dev) <- c(M, N)
    gam_dev <- t(gam_dev)
    gamI <- SqrtMeanInverse(t(gam))
    gamI_dev <- gradient(gamI, 1 / (M - 1))

    mq[, r + 1] <- stats::approx(
      time, mq[, r],
      xout = (time[M] - time[1]) * gamI + time[1]
    )$y * sqrt(gamI_dev)

    for (n in 1:N) {
      q[, n, r + 1] <- stats::approx(
        time, q[, n, r],
        xout = (time[M] - time[1]) * gamI + time[1]
      )$y * sqrt(gamI_dev)
      f[, n, r + 1] <- stats::approx(
        time, f[, n, r],
        xout = (time[M] - time[1]) * gamI + time[1]
      )$y
      gam[n, ] <- stats::approx(
        time, gam[n, ],
        xout = (time[M] - time[1]) * gamI + time[1]
      )$y
    }

    qun[r + 1] <- pvecnorm(mq[, r + 1] - mq[, r], 2) / pvecnorm(mq[, r], 2)
  }

  # Aligned data & stats
  fn <- f[, , r + 1]
  qn <- q[, , r + 1]
  q0 <- q[, , 1]
  mean_f0 <- rowMeans(f[, , 1])
  std_f0 <- apply(f[, , 1], 1, stats::sd)
  mean_fn <- rowMeans(fn)
  std_fn <- apply(fn, 1, stats::sd)
  mqn <- mq[, r + 1]

  if (centroid_type == "mean")
    fmean = mean(f0[1, ]) + as.numeric(cumtrapz(time, mqn * abs(mqn)))
  else
    fmean = stats::median(f0[1, ]) + as.numeric(cumtrapz(time, mqn * abs(mqn)))

  gam <- t(gam)
  fgam <- matrix(0, M, N)

  for (n in 1:N)
    fgam[, n] <- stats::approx(
      time, fmean,
      xout = (time[M] - time[1]) * gam[, n] + time[1]
    )$y

  var_fgam <- apply(fgam, 1, stats::var)
  orig_var <- trapz(time, std_f0 ^ 2)
  amp_var <- trapz(time, std_fn ^ 2)
  phase_var <- trapz(time, var_fgam)

  out <- list(
    time = time,
    f0 = f[, , 1],
    q0 = q0,
    fn = fn,
    qn = qn,
    fmean = fmean,
    mqn = mqn,
    warping_functions = gam,
    original_variance = orig_var,
    amplitude_variance = amp_var,
    phase_variance = phase_var,
    qun = qun[1:r],
    inverse_average_warping_function = gamI,
    rsamps = FALSE,
    call = list(
      lambda = lambda,
      penalty_method = penalty_method,
      centroid_type = centroid_type,
      center_warpings = center_warpings,
      smooth_data = smooth_data,
      sparam = sparam,
      parallel = parallel,
      optim_method = optim_method,
      max_iter = max_iter
    )
  )

  class(out) <- "fdawarp"

  if (parallel) parallel::stopCluster(cl)

  out
}
