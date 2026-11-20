#' Elastic Functional Statistical Process Control
#'
#' Fits PCA-based process monitoring models (Hotelling \eqn{T^2} and SPE/Q
#' charts) to aligned functional data in the SRVF framework, with separate
#' models for amplitude (vertical fPCA, [vertFPCA()]) and phase (horizontal
#' fPCA, [horizFPCA()]). The monitoring statistics and control limits follow
#' Ge and Song (2013), Sec. 2.1, Eqs. (2.3)-(2.8), applied to the elastic
#' fPCA features. New functions are monitored with [spc_monitor()].
#'
#' @param warp_data An object of class `fdawarp` from [time_warping()] holding
#'   the aligned normal-operating-condition data.
#' @param var_exp Cumulative proportion of variance used to choose the number
#'   of components: the smallest number whose cumulative proportion reaches
#'   `var_exp`. A scalar, or a vector `c(amplitude, phase)`. Defaults to
#'   `0.95`.
#' @param no Number of components, a scalar or `c(amplitude, phase)`. Where
#'   given (not `NA`), it overrides `var_exp`. Defaults to `NA`.
#' @param alpha Significance level of the control limits. Defaults to `0.01`.
#' @param t2limit Formula of the \eqn{T^2} limit: `"standard"`,
#'   \eqn{k(n^2-1)/(n(n-k)) F}, or `"book"`, \eqn{k(n-1)/(n-k) F}. Defaults to
#'   `"standard"`.
#' @param id Index of the time point used for \eqn{f(id)} in the amplitude
#'   model. Defaults to the midpoint.
#' @param realign A Boolean specifying whether to compute the training
#'   \eqn{T^2}/SPE and the limits from the training functions re-aligned to
#'   the template with [align_new()], the same path [spc_monitor()] uses.
#'   Re-alignment perturbs the SRVFs, so limits fit to the in-sample
#'   residuals are too tight for monitored data; `FALSE` uses the in-sample
#'   features. Defaults to `TRUE`.
#'
#' @return An object of class `fdaspc`, a list with
#' \item{warp_data}{the `fdawarp` object}
#' \item{vpca}{the amplitude model, a `vfpca` object}
#' \item{hpca}{the phase model, an `hfpca` object}
#' \item{alpha, t2limit}{the significance level and \eqn{T^2} limit formula}
#' \item{amp, ph}{the amplitude and phase models: number of components `k`,
#'   number of training functions `n`, limits `T2lim` and `SPElim`, the
#'   chi-square scale `g` and degrees of freedom `h` of the SPE limit, and
#'   the training statistics `T2` and `SPE`}
#' \item{mon}{the results of the last [spc_monitor()] call (`NULL` until
#'   then)}
#'
#' @references Tucker, J. D. (2016). Functional statistical process control
#'   using elastic methods. Proceedings of the Joint Statistical Meetings.
#' @references Ge, Z., Song, Z. (2013). Multivariate Statistical Process
#'   Control. Springer.
#'
#' @export
#' @examples
#' \dontrun{
#'   warp <- time_warping(simu_data$f, simu_data$time)
#'   spc <- fdaspc(warp, alpha = 0.05)
#'   spc <- spc_monitor(spc, 1.5 * simu_data$f[, 1:3])
#'   spc$mon$alarm
#'   plot(spc)
#'   plot_contrib(spc, 1)
#' }
fdaspc <- function(warp_data,
                   var_exp = 0.95,
                   no = NA,
                   alpha = 0.01,
                   t2limit = c("standard", "book"),
                   id = round(length(warp_data$time) / 2),
                   realign = TRUE) {
  if (!inherits(warp_data, "fdawarp") || is.null(warp_data$fn))
    stop("Please align the data with time_warping() first")
  t2limit <- rlang::arg_match(t2limit)
  if (!(alpha > 0 && alpha < 1))
    stop("alpha must be in (0, 1)")
  var_exp <- spc_pair(var_exp)
  no <- spc_pair(no)

  # amplitude: vertical fPCA on [qn; sign(f(id)) sqrt|f(id)|]
  Xin <- spc_amp_features(warp_data$fn, warp_data$qn, id)
  no[1] <- spc_ncomp(Xin, var_exp[1], no[1])
  vpca <- vertFPCA(warp_data, no = no[1], id = id, showplot = FALSE)

  # phase: horizontal fPCA on the shooting vectors
  vec <- SqrtMean(warp_data$warping_functions)$vec
  no[2] <- spc_ncomp(vec, var_exp[2], no[2])
  hpca <- horizFPCA(warp_data, no = no[2], showplot = FALSE)

  # features the limits are fit to
  if (realign) {
    al <- align_new(warp_data, warp_data$f0)
    X <- spc_amp_features(al$fn, al$qn, id)
    V <- spc_ph_features(hpca$mu, al$gam)
  } else {
    X <- Xin
    V <- hpca$vec
  }

  out <- list(
    warp_data = warp_data,
    vpca = vpca,
    hpca = hpca,
    alpha = alpha,
    t2limit = t2limit,
    amp = spc_fit_limits(X, vpca$mqn, vpca$U, vpca$latent, alpha, t2limit),
    ph = spc_fit_limits(V, hpca$vm, hpca$U, hpca$latent, alpha, t2limit),
    mon = NULL
  )
  class(out) <- "fdaspc"
  out
}

#' Monitor New Functions with an Elastic SPC Model
#'
#' Aligns new functions to the training template of an [fdaspc()] model (same
#' lambda, method and penalty) and computes the scores, \eqn{T^2}, SPE,
#' alarms and per-time-point contributions for amplitude and phase. A fault
#' is detected when \eqn{T^2} exceeds its limit or SPE exceeds its limit.
#'
#' @param x An object of class `fdaspc`.
#' @param f_new A matrix (\eqn{M \times N}) of \eqn{N} new functions sampled
#'   at `x$warp_data$time`.
#'
#' @return The `fdaspc` object with `mon` filled in, a list with
#' \item{fn, qn, gam}{the aligned functions, SRVFs and warping functions}
#' \item{amp, ph}{lists with the scores `T`, `T2`, `SPE`, the alarms
#'   `alarmT2`, `alarmSPE` and `alarm`, the limits `T2lim` and `SPElim`, and
#'   the contributions `contT2` and `contSPE` (rows are the new functions and
#'   columns the time points; the amplitude ones have the \eqn{f(id)} term
#'   last)}
#' \item{alarm}{a logical vector, any amplitude or phase alarm}
#'
#' @export
spc_monitor <- function(x, f_new) {
  if (!inherits(x, "fdaspc"))
    stop("`x` must be an fdaspc object")
  al <- align_new(x$warp_data, f_new)

  X <- spc_amp_features(al$fn, al$qn, x$vpca$id)
  amp <- spc_monitor_stats(X, x$vpca$mqn, x$vpca$U, x$vpca$latent, x$amp)
  V <- spc_ph_features(x$hpca$mu, al$gam)
  ph <- spc_monitor_stats(V, x$hpca$vm, x$hpca$U, x$hpca$latent, x$ph)

  x$mon <- list(fn = al$fn, qn = al$qn, gam = al$gam, amp = amp, ph = ph,
                alarm = amp$alarm | ph$alarm)
  x
}

#' Control Charts of an Elastic SPC Model
#'
#' Plots the \eqn{T^2} and SPE charts for amplitude and phase: the training
#' data, followed by the data of the last [spc_monitor()] call. Dashed lines
#' are the control limits and red points are alarms.
#'
#' @param x An object of class `fdaspc`.
#' @param ... Unused.
#'
#' @return Invisibly, `x`.
#' @export
plot.fdaspc <- function(x, ...) {
  mods <- list(x$amp, x$ph)
  news <- if (is.null(x$mon)) list(NULL, NULL) else list(x$mon$amp, x$mon$ph)
  names_m <- c("Amplitude", "Phase")
  stats_s <- c("T2", "SPE")

  oldpar <- graphics::par(mfrow = c(2, 2))
  on.exit(graphics::par(oldpar))
  for (a in 1:2) {
    for (s in 1:2) {
      tr <- mods[[a]][[stats_s[s]]]
      lim <- mods[[a]][[paste0(stats_s[s], "lim")]]
      ntr <- length(tr)
      nw <- if (is.null(news[[a]])) numeric(0) else news[[a]][[stats_s[s]]]
      idx <- ntr + seq_along(nw)
      ylab <- if (s == 1) expression(T^2) else "SPE"
      main <- if (s == 1) {
        bquote(.(names_m[a]) ~ T^2 ~ "(k =" ~ .(mods[[a]]$k) * ")")
      } else {
        sprintf("%s SPE (k = %d)", names_m[a], mods[[a]]$k)
      }
      plot(seq_len(ntr), tr, pch = 20, col = "grey50",
           xlim = c(1, ntr + length(nw)), ylim = range(0, tr, nw, lim),
           xlab = "Observation", ylab = ylab, main = main)
      if (length(nw) > 0) {
        al <- nw > lim
        graphics::lines(idx, nw, type = "o", pch = 20, col = "blue")
        graphics::points(idx[al], nw[al], pch = 21, col = "red", bg = "red")
        graphics::abline(v = ntr + 0.5, lty = 3)
      }
      graphics::abline(h = lim, lty = 2, lwd = 1.5)
    }
  }
  invisible(x)
}

#' Contributions of a Monitored Function
#'
#' Plots the \eqn{T^2} and SPE contributions of each time point for the
#' amplitude and phase models, to locate where in time a fault occurs.
#'
#' @param x An object of class `fdaspc` after [spc_monitor()].
#' @param idx Index of the function passed to the last [spc_monitor()] call.
#'
#' @return Invisibly, `x`.
#' @export
plot_contrib <- function(x, idx) {
  if (!inherits(x, "fdaspc") || is.null(x$mon))
    stop("Please call spc_monitor() first")
  t <- x$warp_data$time
  M <- length(t)
  oldpar <- graphics::par(mfrow = c(2, 2))
  on.exit(graphics::par(oldpar))
  plot(t, x$mon$amp$contT2[idx, 1:M], type = "l", lwd = 1.5, xlab = "t",
       ylab = "", main = expression("Amplitude " * T^2 * " contribution"))
  plot(t, x$mon$amp$contSPE[idx, 1:M], type = "l", lwd = 1.5, xlab = "t",
       ylab = "", main = "Amplitude SPE contribution")
  plot(t, x$mon$ph$contT2[idx, ], type = "l", lwd = 1.5, xlab = "t",
       ylab = "", main = expression("Phase " * T^2 * " contribution"))
  plot(t, x$mon$ph$contSPE[idx, ], type = "l", lwd = 1.5, xlab = "t",
       ylab = "", main = "Phase SPE contribution")
  invisible(x)
}

# scalar -> c(x, x), keep a 2-vector
spc_pair <- function(x) {
  if (length(x) == 1) return(c(x, x))
  if (length(x) == 2) return(as.numeric(x))
  stop("expected a scalar or a 2-element c(amplitude, phase) vector")
}

# number of components: `no` when given, else the first whose cumulative
# proportion of variance reaches var_exp
spc_ncomp <- function(X, var_exp, no) {
  if (!is.na(no)) return(as.integer(no))
  if (var_exp > 1) stop("var_exp is greater than 1")
  s <- svd(stats::cov(t(X)), nu = 0, nv = 0)$d
  which(cumsum(s) / sum(s) >= var_exp)[1]
}

# columns [qn; sign(f(id)) sqrt|f(id)|] as in vertFPCA()
spc_amp_features <- function(fn, qn, id) {
  m <- sign(fn[id, ]) * sqrt(abs(fn[id, ]))
  rbind(qn, m)
}

# shooting vectors at the Karcher mean, the same map as SqrtMean()
spc_ph_features <- function(mu, gam) {
  M <- nrow(gam)
  binsize <- 1 / (M - 1)
  apply(gam, 2, function(g)
    inv_exp_map(mu, sqrt(pmax(gradient(g, binsize), 0))))
}

# Ge & Song Eqs. (2.3)-(2.6) and contributions; X is p x N
spc_scores <- function(X, mu, U, latent) {
  U <- matrix(U, nrow = nrow(X))
  Xc <- t(X - mu)
  Tsc <- Xc %*% U
  E <- Xc - Tsc %*% t(U)
  T2 <- rowSums(sweep(Tsc ^ 2, 2, latent, "/"))
  SPE <- rowSums(E ^ 2)
  D <- sweep(U, 2, sqrt(latent), "/")
  contT2 <- (Xc %*% D) %*% t(D) * Xc # rows sum to T2
  contSPE <- E ^ 2                   # rows sum to SPE
  list(T = Tsc, T2 = T2, SPE = SPE, contT2 = contT2, contSPE = contSPE)
}

# control limits from the training data, Ge & Song Eqs. (2.7)-(2.8)
spc_fit_limits <- function(X, mu, U, latent, alpha, t2limit) {
  n <- ncol(X)
  k <- length(latent)
  if (k >= n)
    stop(sprintf(paste("Need fewer components than functions for the T^2",
                       "limit (k = %d, n = %d)."), k, n))
  S <- spc_scores(X, mu, U, latent)

  Fcrit <- stats::qf(1 - alpha, k, n - k)
  T2lim <- if (t2limit == "book") {
    k * (n - 1) / (n - k) * Fcrit
  } else {
    k * (n ^ 2 - 1) / (n * (n - k)) * Fcrit
  }

  mSPE <- mean(S$SPE)
  vSPE <- stats::var(S$SPE)
  if (vSPE < .Machine$double.eps) {
    g <- .Machine$double.eps
    h <- 1
  } else {
    g <- vSPE / (2 * mSPE)
    h <- 2 * mSPE ^ 2 / vSPE
  }
  SPElim <- g * stats::qchisq(1 - alpha, h)

  list(k = k, n = n, T2lim = T2lim, SPElim = SPElim, g = g, h = h,
       T2 = S$T2, SPE = S$SPE)
}

spc_monitor_stats <- function(X, mu, U, latent, mdl) {
  out <- spc_scores(X, mu, U, latent)
  out$alarmT2 <- out$T2 > mdl$T2lim
  out$alarmSPE <- out$SPE > mdl$SPElim
  out$alarm <- out$alarmT2 | out$alarmSPE
  out$T2lim <- mdl$T2lim
  out$SPElim <- mdl$SPElim
  out
}
