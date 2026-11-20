#' @importFrom graphics plot
#' @export
plot.fdawarp <- function(x, ...){

  colors = c(
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#B3B3B3"
  )

  M <- nrow(x$fn)
  mean_f0 <- rowMeans(x$f0)
  std_f0 <- apply(x$f0, 1, stats::sd)
  mean_fn <- rowMeans(x$fn)
  std_fn <- apply(x$fn, 1, stats::sd)

  graphics::matplot(
    x = x$time,
    y = x$f0,
    type = "l",
    col = colors,
    main = "Original Data",
    xlab ="Time",
    ylab = "Amplitude"
  )

  graphics::matplot(
    x = (0:(M - 1)) / (M - 1),
    y = x$warping_functions,
    type = "l",
    col = colors,
    main = "Warping functions",
    xlab = "Time",
    ylab = "Warped Time"
  )

  graphics::matplot(
    x = x$time,
    y = x$fn,
    type = "l",
    col = colors,
    main = bquote(paste("Warped Data (", lambda == .(x$call$lambda), ")")),
    xlab ="Time",
    ylab = "Amplitude"
  )

  graphics::matplot(
    x = x$time,
    y = cbind(mean_f0, mean_f0 + std_f0, mean_f0 - std_f0),
    type = "l",
    lty = 1,
    col = c("black", "#66C2A5", "#66C2A5"),
    main = bquote(paste("Original Data: ", Mean %+-% STD)),
    xlab ="Time",
    ylab = "Amplitude"
  )

  graphics::matplot(
    x = x$time,
    y = cbind(mean_fn, mean_fn + std_fn, mean_fn - std_fn),
    type = "l",
    lty = 1,
    col = c("black", "#66C2A5", "#66C2A5"),
    main = bquote(paste(
      "Warped Data: ", Mean %+-% STD,
      " (", lambda == .(x$call$lambda), ")"
    )),
    xlab ="Time",
    ylab = "Amplitude"
  )

  if (x$call$centroid_type == "mean")
    plot(
      x = x$time,
      y = x$fmean,
      type = "l",
      col = "#66C2A5",
      main = bquote(paste(f[mean], " (", lambda == .(x$call$lambda), ")")),
      xlab ="Time"
      ,
      ylab = "Amplitude"
    )
  else
    plot(
      x = x$time,
      y = x$fmean,
      type = "l",
      col = "#66C2A5",
      main = bquote(paste(f[median], " (", lambda == .(x$call$lambda), ")")),
      xlab = "Time",
      ylab = "Amplitude"
    )

  # MCMC diagnostics when the Bayesian alignment was used
  if (is.list(x$mcmc) && (!is.null(x$mcmc$sigma1) || !is.null(x$mcmc$sigma)))
    plot_mcmc(x)
}

#' Plot MCMC Diagnostics of a Bayesian Alignment
#'
#' Plots trace, posterior and mixing diagnostics of the chain stored in the
#' `mcmc` component of an `fdawarp` object returned by [time_warping_bayes()]
#' (pCN) or [time_warping_bayes_infHMC()]. Plots that need the extra
#' information (`extrainfo = TRUE`) are drawn only when it is available.
#'
#' The plots are:
#' - chain diagnostics: trace of \eqn{\sigma} and of the log-likelihood with
#'   the burn-in marked, posterior of \eqn{\sigma} with its 95% interval,
#'   autocorrelation of \eqn{\sigma} with the effective sample size, running
#'   acceptance rate, and proposals and acceptance by pCN step (pCN) or the
#'   HMC step size of each function (infHMC);
#' - trace of the norm of each function's warping coefficients;
#' - posterior template with a 95% credible band;
#' - posterior warping functions with 95% credible bands.
#'
#' @param x An object of class `fdawarp` with an `mcmc` component.
#'
#' @return Invisibly, `x`.
#' @export
#' @examples
#' \dontrun{
#'   out <- time_warping_bayes_infHMC(simu_data$f, simu_data$time)
#'   plot_mcmc(out)
#' }
plot_mcmc <- function(x) {
  mc <- x$mcmc
  if (!is.list(mc) || (is.null(mc$sigma1) && is.null(mc$sigma)))
    stop("No MCMC output available, run time_warping_bayes() or ",
         "time_warping_bayes_infHMC() first")
  is_hmc <- !is.null(mc$epsilon)

  sigma1 <- as.numeric(if (!is.null(mc$sigma1)) mc$sigma1 else mc$sigma)
  niter <- length(sigma1)
  burnin <- if (!is.null(mc$burnin)) {
    mc$burnin
  } else if (!is.null(x$call$burnin)) {
    x$call$burnin
  } else if (!is.null(mc$q_star_mat)) {
    niter - ncol(mc$q_star_mat) + 1
  } else {
    floor(niter / 2)
  }
  burnin <- min(max(burnin, 1), niter)
  logl <- as.numeric(if (!is.null(mc$logl)) mc$logl else x$qun)
  has_extra <- length(mc$accept) > 0
  post <- burnin:niter
  col <- "#66C2A5"
  mark_burnin <- function() graphics::abline(v = burnin, lty = 2)
  not_stored <- function(what, main) {
    graphics::plot.new()
    graphics::title(main = main)
    graphics::text(0.5, 0.5, paste0(what, " not stored\n(set extrainfo = TRUE)"))
  }

  # chain diagnostics -------------------------------------------------------
  oldpar <- graphics::par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))
  on.exit(graphics::par(oldpar), add = TRUE)

  plot(1:niter, sigma1, type = "l", col = col, xlab = "iteration",
       ylab = expression(sigma), main = expression("Trace of " * sigma))
  mark_burnin()

  if (length(logl) == niter) {
    plot(1:niter, logl, type = "l", col = col, xlab = "iteration",
         ylab = "log-likelihood", main = "Trace of log-likelihood")
    mark_burnin()
  } else {
    not_stored("log-likelihood", "Trace of log-likelihood")
  }

  sp <- sigma1[post]
  graphics::hist(sp, freq = FALSE, col = col, border = NA,
                 xlab = expression(sigma), ylab = "density",
                 main = expression("Posterior of " * sigma * " (95% CI)"))
  graphics::abline(v = mean(sp), lwd = 1.5)
  graphics::abline(v = stats::quantile(sp, c(0.025, 0.975)), lty = 2)

  maxlag <- min(50, length(post) - 1)
  rho <- mcmc_acf(sp, maxlag)
  ess <- length(post) / max(1 + 2 * sum(rho[-1]), 1)
  plot(0:maxlag, rho, type = "h", col = col, xlab = "lag", ylab = "ACF",
       main = bquote("ACF of " * sigma * " (ESS" %~~% .(round(ess)) * ")"))
  graphics::abline(h = c(-1, 1) * 1.96 / sqrt(length(post)), lty = 2)
  graphics::abline(h = 0)

  if (has_extra) {
    # HMC stores (N, iter - 1) indicators: average over functions
    acc <- if (is_hmc) colMeans(matrix(as.numeric(mc$accept),
                                       nrow = NROW(mc$accept)))
           else as.numeric(mc$accept)
    lab <- if (is_hmc) "warps" else "q*"
    it <- seq_along(acc) + 1
    plot(it, cumsum(acc) / seq_along(acc), type = "l", col = col,
         ylim = c(0, 1), xlab = "iteration", ylab = "rate",
         main = sprintf("Running acceptance of %s (%.2f)", lab, mean(acc)))
    mark_burnin()
  } else {
    not_stored("acceptance", "Running acceptance of q*")
  }

  if (is_hmc) {
    eps <- t(mc$epsilon)
    graphics::matplot(1:nrow(eps), eps, type = "l", lty = 1, log = "y",
                      xlab = "iteration", ylab = expression(epsilon),
                      main = "HMC step size per function")
    mark_burnin()
  } else if (has_extra && !is.null(mc$betas_ind)) {
    bind <- as.integer(mc$betas_ind)
    acc <- as.logical(mc$accept)
    nb <- if (!is.null(mc$zpcn)) length(mc$zpcn$betas) else max(bind)
    ok <- bind > 0
    nprop <- tabulate(bind[ok], nbins = nb)
    nacc <- tabulate(bind[ok & acc], nbins = nb)
    names_arg <- if (!is.null(mc$zpcn)) format(mc$zpcn$betas) else 1:nb
    bp <- graphics::barplot(rbind(nprop - nacc, nacc), names.arg = names_arg,
                            col = c("#FC8D62", col), border = NA,
                            ylim = c(0, 1.15 * max(c(nprop, 1))),
                            xlab = if (!is.null(mc$zpcn)) "pCN beta"
                                   else "pCN step index",
                            ylab = "proposals",
                            main = "q* proposals by pCN step (acc. rate)",
                            legend.text = c("rejected", "accepted"),
                            args.legend = list(x = "topright", bty = "n"))
    graphics::text(bp, nprop, sprintf("%.2f", nacc / pmax(nprop, 1)), pos = 3)
  } else {
    not_stored("pCN step usage", "pCN step usage")
  }
  graphics::mtext("MCMC diagnostics", outer = TRUE, cex = 1.2)
  graphics::par(oldpar)

  # warping coefficient traces ---------------------------------------------
  if (length(mc$g_coef) > 0) {
    # pCN stores (iter, D, N)
    gnorm <- sqrt(apply(mc$g_coef ^ 2, c(1, 3), sum))
    graphics::matplot(1:nrow(gnorm), gnorm, type = "l", lty = 1,
                      xlab = "iteration",
                      ylab = "norm of warping coefficients",
                      main = "Trace of warping coefficient norm per function")
    mark_burnin()
  } else if (length(mc$v_coef) > 0) {
    # infHMC stores (D, N, nsamp) post burn-in draws only
    vnorm <- sqrt(apply(mc$v_coef ^ 2, c(3, 2), sum))
    graphics::matplot(1:nrow(vnorm), vnorm, type = "l", lty = 1,
                      xlab = "posterior draw (after burn-in)",
                      ylab = "norm of warping coefficients",
                      main = "Trace of warping coefficient norm per function")
  }

  if (is.null(mc$fmean_mat) || is.null(mc$gamma_mat))
    return(invisible(x))

  # posterior template ------------------------------------------------------
  M <- length(x$time)
  ts <- seq(x$time[1], x$time[M], length.out = nrow(mc$fmean_mat))
  ndraw <- ncol(mc$fmean_mat)
  idx <- unique(round(seq(1, ndraw, length.out = min(ndraw, 100))))
  graphics::matplot(ts, mc$fmean_mat[, idx, drop = FALSE], type = "l",
                    lty = 1, col = grDevices::adjustcolor("grey70", 0.3),
                    xlab = "Time", ylab = "Amplitude",
                    ylim = range(mc$fmean_mat, x$fmean),
                    main = "Posterior template with 95% credible band")
  graphics::polygon(c(ts, rev(ts)),
                    c(mc$fmean_stats[, 1], rev(mc$fmean_stats[, 2])),
                    col = grDevices::adjustcolor(col, 0.35), border = NA)
  graphics::lines(x$time, x$fmean, lwd = 2)
  graphics::legend("topright", "posterior mean", lwd = 2, bty = "n")

  # posterior warping functions --------------------------------------------
  N <- dim(mc$gamma_mat)[3]
  tg <- seq(0, 1, length.out = dim(mc$gamma_mat)[1])
  cols <- grDevices::hcl.colors(N, "Dark 3")
  plot(c(0, 1), c(0, 1), type = "l", lty = 3, asp = 1, xlab = "Time",
       ylab = "Warped Time",
       main = "Posterior warping functions with 95% credible bands")
  for (jj in 1:N) {
    band <- t(apply(mc$gamma_mat[, , jj, drop = FALSE], 1, stats::quantile,
                    probs = c(0.025, 0.975)))
    graphics::polygon(c(tg, rev(tg)), c(band[, 1], rev(band[, 2])),
                      col = grDevices::adjustcolor(cols[jj], 0.2), border = NA)
  }
  tw <- seq(0, 1, length.out = nrow(x$warping_functions))
  for (jj in 1:N)
    graphics::lines(tw, x$warping_functions[, jj], col = cols[jj])

  invisible(x)
}

# sample autocorrelation of x for lags 0:maxlag
mcmc_acf <- function(x, maxlag) {
  x <- as.numeric(x) - mean(x)
  n <- length(x)
  c0 <- sum(x ^ 2)
  rho <- numeric(maxlag + 1)
  if (c0 == 0) {
    rho[1] <- 1
    return(rho)
  }
  for (k in 0:maxlag)
    rho[k + 1] <- sum(x[1:(n - k)] * x[(1 + k):n]) / c0
  rho
}
