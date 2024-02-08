#' Functional Boxplot
#'
#' This function computes the required statistics for building up a boxplot of
#' the aligned functional data. Since the process of alignment provides
#' separation of phase and amplitude variability, the computed boxplot can focus
#' either on amplitude variability or phase variability.
#'
#' The function [boxplot.fdawarp()] returns optionally an object of class either
#' `ampbox` if `variability_type = "amplitude"` or `phbox` if `variability_type
#' = "phase"`. `S3` methods specialized for objects of these classes are
#' provided as well to avoid re-computation of the boxplot statistics.
#'
#' @param x An object of class `fdawarp` typically produced by [time_warping()]
#'   or of class `ampbox` or `phbox` typically produced by [boxplot.fdawarp()].
#' @param variability_type A string specifying which kind of variability should
#'   be displayed in the boxplot. Choices are `"amplitude"` or `"phase"`.
#'   Defaults to `"amplitude"`.
#' @param alpha A numeric value specifying the quantile value. Defaults to
#'   \eqn{0.05} which uses the \eqn{95\%} quantile.
#' @param range A positive numeric value specifying how far the plot whiskers
#'   extend out from the box. The whiskers extend to the most extreme data point
#'   which is no more than `range` times the interquartile range from the box.
#'   Defaults to `1.0`.
#' @param what A string specifying what the function should return. Choices are
#'   `"plot"`, `"stats"` or `"plot+stats"`. Defaults to `"stats"`.
#' @param ... Unused here.
#'
#' @return If `what` contains `stats`, a list containing the computed statistics
#'   necessary for drawing the boxplot. Otherwise, the function simply draws the
#'   boxplot and no object is returned.
#'
#' @importFrom graphics boxplot
#' @export
#'
#' @examples
#' \dontrun{
#' out <- time_warping(simu_data$f, simu_data$time)
#' boxplot(out, what = "stats")
#' }
boxplot.fdawarp <- function(x,
                            variability_type = c("amplitude", "phase"),
                            alpha = 0.05,
                            range = 1.0,
                            what = c("stats", "plot", "plot+stats"),
                            ...) {
  variability_type <- rlang::arg_match(variability_type)
  what <- rlang::arg_match(what)

  if (x$call$centroid_type != "median") {
    cli::cli_alert_warning(
      "The argument {.arg x} is of class {.cls fdawarp} but has not been
      computed using the median as centroid type."
    )
    cli::cli_alert_info(
      'Rerunning {.fn time_warping} with {.code centroid_type = "median"}...'
    )
    x <- time_warping(
      f = x$f0, time = x$time,
      lambda = x$call$lambda,
      penalty_method = x$call$penalty_method,
      centroid_type = "median",
      center_warpings = x$call$center_warpings,
      smooth_data = x$call$smooth_data,
      sparam = x$call$sparam,
      parallel = x$call$parallel,
      optim_method = x$call$optim_method,
      max_iter = x$call$max_iter
    )
  }

  plot_data <- if (variability_type == "amplitude")
    ampbox_data(x, alpha = alpha, ka = range)
  else
    phbox_data(x, alpha = alpha, kp = range)

  if (what == "plot" || what == "plot+stats")
    boxplot(plot_data)

  if (what == "stats" || what == "plot+stats")
    return(plot_data)
}

#' @rdname boxplot.fdawarp
#' @export
boxplot.ampbox <- function(x, ...) {
  fmedian <- x$fmedian
  maxx <- x$maxx
  minn <- x$minn
  Q1 <- x$Q1
  Q1a <- x$Q1a
  Q3 <- x$Q3
  Q3a <- x$Q3a
  time <- x$time
  M <- length(fmedian)
  ymin <- min(c(min(fmedian), min(Q1), min(Q3), min(maxx), min(minn)))
  ymax <- max(c(max(fmedian), max(Q1), max(Q3), max(maxx), max(minn)))
  plot(
    time, fmedian,
    col = "black",
    xlab = "Time",
    main = "Amplitude Boxplot",
    type = "l",
    ylim = c(ymin, ymax)
  )
  graphics::lines(time, Q1, col = "blue")
  graphics::lines(time, Q3, col = "blue")
  graphics::lines(time, Q1a, col = "green")
  graphics::lines(time, Q3a, col = "green")
  graphics::lines(time, maxx, col = "red")
  graphics::lines(time, minn, col = "red")

  s <- seq(0, 1, length.out = 100)
  Fs2 <- matrix(0, length(time), 595)
  Fs2[, 1] <- (1 - s[1]) * minn + s[1] * Q1
  for (j in 2:100) {
    Fs2[, j] <- (1 - s[j]) * minn + s[j] * Q1a
    Fs2[, 99 + j] <- (1 - s[j]) * Q1a + s[j] * Q1
    Fs2[, 198 + j] <- (1 - s[j]) * Q1 + s[j] * fmedian
    Fs2[, 297 + j] <- (1 - s[j]) * fmedian + s[j] * Q3
    Fs2[, 396 + j] <- (1 - s[j]) * Q3 + s[j] * Q3a
    Fs2[, 495 + j] <- (1 - s[j]) * Q3a + s[j] * maxx
  }
  d1 <- sqrt(trapz(time, (x$qmedian - x$Q1_q) ^ 2))
  d1a <- sqrt(trapz(time, (x$Q1_q - x$Q1a_q) ^ 2))
  dl <- sqrt(trapz(time, (x$Q1a_q - x$min_q) ^ 2))
  d3 <- sqrt(trapz(time, (x$qmedian - x$Q3_q) ^ 2))
  d3a <- sqrt(trapz(time, (x$Q3_q - x$Q3a_q) ^ 2))
  du <- sqrt(trapz(time, (x$Q3a_q - x$max_q) ^ 2))
  part1 <- seq(-d1 - d1a - dl, -d1 - d1a, length.out = 100)
  part2 <- seq(-d1 - d1a, -d1, length.out = 100)
  part3 <- seq(-d1, 0, length.out = 100)
  part4 <- seq(0, d3, length.out = 100)
  part5 <- seq(d3, d3 + d3a, length.out = 100)
  part6 <- seq(d3 + d3a, d3 + d3a + du, length.out = 100)
  allparts <- c(
    part1, part2[2:100], part3[2:100],
    part4[2:100], part5[2:100], part6[2:100]
  )

  if (requireNamespace("plot3Drgl", quietly = TRUE)) {
    plot3D::persp3D(
      x = time,
      y = allparts,
      z = Fs2,
      col = viridisLite::viridis(128),
      plot = FALSE,
      main = "Amplitude Surface Plot",
      ticktype = "detailed",
      box = FALSE
    ) +
      plot3D::lines3D(
        x = time,
        y = rep(0, M),
        z = fmedian,
        col = "black",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(-d1, M),
        z = Q1,
        col = "blue",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(-d1 - d1a, M),
        z = Q1a,
        col = "green",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(-d1 - d1a - dl, M),
        z = minn,
        col = "red",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(d3, M),
        z = Q3,
        col = "blue",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(d3 + d3a, M),
        z = Q3a,
        col = "green",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(d3 + d3a + du, M),
        z = maxx,
        col = "red",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      )
    plot3Drgl::plotrgl()
    rgl::par3d("windowRect" = c(0, 0, 640, 640))
    rgl::grid3d(c("x", "y+", "z"))
    rgl::axes3d(c('x--', "y--", 'z'))
    rgl::title3d(xlab = "Time", ylab = "Distance")
  } else {
    graphics::image(
      time,
      allparts,
      Fs2,
      main = "Surface Plot",
      ylab = "",
      col = viridisLite::viridis(128)
    )
    graphics::lines(time, rep(0, M), col = "black", lwd = 1)
    graphics::lines(time, rep(-d1, M), col = "blue", lwd = 1)
    graphics::lines(time, rep(-d1 - d1a, M), col = "green", lwd = 1)
    graphics::lines(time, rep(-d1 - d1a - dl, M), col = "red", lwd = 1)
    graphics::lines(time, rep(d3, M), col = "blue", lwd = 1)
    graphics::lines(time, rep(d3 + d3a, M), col = "green", lwd = 1)
    graphics::lines(time, rep(d3 + d3a + du, M), col = "red", lwd = 1)
  }
}

#' Amplitude Boxplot Data
#'
#' This function constructs the amplitude boxplot.
#'
#' @param warp_median fdawarp object from [time_warping()] of aligned data using
#'   the median
#' @param alpha quantile value (default=.05, i.e., 95%)
#' @param ka scalar for outlier cutoff (default=1)
#'
#' @return Returns an `ampbox` object containing:
#'
#' - `median_y`: median function
#' - `Q1`: First quartile
#' - `Q3`: Second quartile
#' - `Q1a`: First quantile based on alpha
#' - `Q3a`: Second quantile based on alpha
#' - `minn`: minimum extreme function
#' - `maxx`: maximum extreme function
#' - `outlier_index`: indexes of outlier functions
#' - `fmedian`: median function
#'
#' @keywords srvf alignment boxplot
#'
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A geometric
#'   approach to visualization of variability in functional data." Journal of
#'   the American Statistical Association in press: 1-34.
#'
#' @keywords internal
ampbox_data <- function(warp_median, alpha = 0.05, ka = 1) {
  fn <- warp_median$fn
  median_y <- warp_median$fmean
  qn <- warp_median$qn
  qmedian <- warp_median$mqn
  time <- warp_median$time

  if (warp_median$rsamps) {
    fn <- warp_median$fs
    qn <- warp_median$qs
  }

  M <- nrow(fn)
  N <- ncol(fn)
  lambda <- 0.5

  # translation
  translation <- rep(0, N)
  for (i in 1:N)
    translation[i] <- trapz(time, fn[, i] / (time[M] - time[1]))

  # compute amplitude distances
  dy <- rep(0, N)
  for (i in 1:N)
    dy[i] <- sqrt(trapz(time, (qmedian - qn[, i])^2))
  dy_ordering <- sort(dy, index.return = TRUE)$ix
  CR_50 <- dy_ordering[1:round(N / 2)]  # 50% central region
  m <- max(dy[CR_50])  # maximal amplitude distance with 50% central region

  # identify amplitude quartiles
  angle <- matrix(0, nrow = length(CR_50), ncol = length(CR_50))
  energy <- matrix(0, nrow = length(CR_50), ncol = length(CR_50))

  for (i in 1:(length(CR_50) - 1)) {
    for (j in (i + 1):length(CR_50)) {
      q1 <- qn[, CR_50[i]] - qmedian
      q3 <- qn[, CR_50[j]] - qmedian
      q1 <- q1 / sqrt(trapz(time, q1 * q1))
      q3 <- q3 / sqrt(trapz(time, q3 * q3))
      angle[i, j] <- trapz(time, q1 * q3)
      energy[i, j] <- (1 - lambda) * (dy[CR_50[i]] / m + dy[CR_50[j]] / m) -
        lambda * (angle[i, j] + 1)
    }
  }

  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1_index <- CR_50[maxloc[1, 1]]
  Q3_index <- CR_50[maxloc[1, 2]]
  Q1_q <- qn[, Q1_index]
  Q3_q <- qn[, Q3_index]
  Q1 <- fn[, Q1_index]
  Q3 <- fn[, Q3_index]

  # identify amplitude quantile
  dy_ordering <- sort(dy, index.return = TRUE)$ix
  CR_alpha <- dy_ordering[1:round(N * (1 - alpha))]  # (1-alpha)% central region
  m <- max(dy[CR_alpha])  # maximal amplitude distance with (1-alpha)% central region
  angle <- matrix(0, nrow = length(CR_alpha), ncol = length(CR_alpha))
  energy <- matrix(0, nrow = length(CR_alpha), ncol = length(CR_alpha))

  for (i in 1:(length(CR_alpha) - 1)) {
    for (j in (i + 1):length(CR_alpha)) {
      q1 <- qn[, CR_alpha[i]] - qmedian
      q3 <- qn[, CR_alpha[j]] - qmedian
      q1 <- q1 / sqrt(trapz(time, q1 * q1))
      q3 <- q3 / sqrt(trapz(time, q3 * q3))
      angle[i, j] <- trapz(time, q1 * q3)
      energy[i, j] <- (1 - lambda) * (dy[CR_alpha[i]] / m + dy[CR_alpha[j]] / m) -
        lambda * (angle[i, j] + 1)
    }
  }

  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1a_index <- CR_alpha[maxloc[1, 1]]
  Q3a_index <- CR_alpha[maxloc[1, 2]]
  Q1a_q <- qn[, Q1a_index]
  Q3a_q <- qn[, Q3a_index]
  Q1a <- fn[, Q1a_index]
  Q3a <- fn[, Q3a_index]

  # compute amplitude whiskers
  IQR <- dy[Q1_index] + dy[Q3_index]
  v1 <- Q1_q - qmedian
  v3 <- Q3_q - qmedian
  upper_q <- Q3_q + ka * IQR * v3 / sqrt(trapz(time, v3 * v3))
  lower_q <- Q1_q + ka * IQR * v1 / sqrt(trapz(time, v1 * v1))
  upper <- cumtrapz(time, upper_q * abs(upper_q))
  lower <- cumtrapz(time, lower_q * abs(lower_q))

  upper_dis <- sqrt(trapz(time, (upper_q - qmedian)^2))
  lower_dis <- sqrt(trapz(time, (lower_q - qmedian)^2))
  whisker_dis <- max(c(upper_dis, lower_dis))

  # identify amplitude outliers
  outlier_index <- c()
  for (i in 1:N) {
    if (dy[dy_ordering[N + 1 - i]] > whisker_dis)
      outlier_index <- c(outlier_index, dy_ordering[N + 1 - i])
    else
      break
  }

  # identify amplitude extremes
  distance_to_upper <- rep(Inf, N)
  distance_to_lower <- rep(Inf, N)
  out_50_CR <- setdiff(setdiff(1:N, CR_50), outlier_index)

  for (i in 1:length(out_50_CR)) {
    j <- out_50_CR[i]
    distance_to_upper[j] = sqrt(trapz(time, (upper_q - qn[, j])^2))
    distance_to_lower[j] = sqrt(trapz(time, (lower_q - qn[, j])^2))
  }

  max_index <- which.min(distance_to_upper)
  min_index <- which.min(distance_to_lower)
  min_q <- qn[, min_index]
  max_q <- qn[, max_index]
  minn <- fn[, min_index]
  maxx <- fn[, max_index]

  out <- list(
    median_y = median_y,
    Q1 = Q1,
    Q3 = Q3,
    Q1a = Q1a,
    Q3a = Q3a,
    minn = minn,
    maxx = maxx,
    outlier_index = outlier_index
  )
  out$fmedian <- median_y
  out$time <- time
  out$qmedian <- qmedian
  out$Q1_q <- Q1_q
  out$Q1a_q <- Q1a_q
  out$Q3_q <- Q3_q
  out$Q3a_q <- Q3a_q
  out$min_q <- min_q
  out$max_q <- max_q
  out$dist <- dy
  out$Q1a_index <- Q1a_index
  out$Q1_index <- Q1_index
  out$Q3a_index <- Q3a_index
  out$Q3_index <- Q3_index

  class(out) <- 'ampbox'

  out
}

#' @rdname boxplot.fdawarp
#' @export
boxplot.phbox <- function(x, ...) {
  median_x <- x$median_x
  maxx <- x$maxx
  minn <- x$minn
  Q1 <- x$Q1
  Q1a <- x$Q1a
  Q3 <- x$Q3
  Q3a <- x$Q3a
  time <- x$time
  M <- length(median_x)
  plot(
    time,
    median_x,
    col = "black",
    xlab = "Time",
    main = "Phase Boxplot",
    type = "l",
    ylim = c(0, 1)
  )
  graphics::lines(time, Q1, col = "blue")
  graphics::lines(time, Q3, col = "blue")
  graphics::lines(time, Q1a, col = "green")
  graphics::lines(time, Q3a, col = "green")
  graphics::lines(time, maxx, col = "red")
  graphics::lines(time, minn, col = "red")

  s <- seq(0, 1, length.out = 100)
  Fs2 <- matrix(0, length(time), 595)
  Fs2[, 1] <- (1 - s[1]) * (minn - time) + s[1] * (Q1 - time)
  for (j in 2:100) {
    Fs2[, j] <- (1 - s[j]) * (minn - time) + s[j] * (Q1a - time)
    Fs2[, 99 + j] <- (1 - s[j]) * (Q1a - time) + s[j] * (Q1 - time)
    Fs2[, 198 + j] <- (1 - s[j]) * (Q1 - time) + s[j] * (median_x - time)
    Fs2[, 297 + j] <- (1 - s[j]) * (median_x - time) + s[j] * (Q3 - time)
    Fs2[, 396 + j] <- (1 - s[j]) * (Q3 - time) + s[j] * (Q3a - time)
    Fs2[, 495 + j] <- (1 - s[j]) * (Q3a - time) + s[j] * (maxx - time)
  }
  d1 <- acos(inner_product(x$psi_median, x$Q1_psi))
  d1a <- acos(inner_product(x$Q1_psi, x$Q1a_psi))
  dl <- acos(inner_product(x$Q1a_psi, x$min_psi))
  d3 <- acos(inner_product(x$psi_median, x$Q3_psi))
  d3a <- acos(inner_product(x$Q3_psi, x$Q3a_psi))
  du <- acos(inner_product(x$Q3a_psi, x$max_psi))
  part1 <- seq(-d1 - d1a - dl, -d1 - d1a, length.out = 100)
  part2 <- seq(-d1 - d1a, -d1, length.out = 100)
  part3 <- seq(-d1, 0, length.out = 100)
  part4 <- seq(0, d3, length.out = 100)
  part5 <- seq(d3, d3 + d3a, length.out = 100)
  part6 <- seq(d3 + d3a, d3 + d3a + du, length.out = 100)
  allparts <- c(
    part1, part2[2:100], part3[2:100],
    part4[2:100], part5[2:100], part6[2:100]
  )

  if (requireNamespace("plot3Drgl", quietly = TRUE)) {
    plot3D::persp3D(
      x = time,
      y = allparts,
      z = Fs2,
      col = viridisLite::viridis(128),
      plot = FALSE,
      main = "Phase Surface Plot",
      ticktype = "detailed",
      box = FALSE
    ) +
      plot3D::lines3D(
        x = time,
        y = rep(0, M),
        z = (median_x - time),
        col = "black",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(-d1, M),
        z = (Q1 - time),
        col = "blue",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(-d1 - d1a, M),
        z = (Q1a - time),
        col = "green",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(-d1 - d1a - dl, M),
        z = (minn - time),
        col = "red",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(d3, M),
        z = (Q3 - time),
        col = "blue",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(d3 + d3a, M),
        z = (Q3a - time),
        col = "green",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      ) +
      plot3D::lines3D(
        x = time,
        y = rep(d3 + d3a + du, M),
        z = (maxx - time),
        col = "red",
        lwd = 6,
        add = TRUE,
        plot = FALSE
      )
    plot3Drgl::plotrgl()
    rgl::par3d("windowRect" = c(0, 0, 640, 640))
    rgl::grid3d(c("x", "y+", "z"))
    rgl::axes3d(c('x--', "y--", 'z'))
    rgl::title3d(xlab = "Time", ylab = "Distance")
  } else {
    graphics::image(
      time,
      allparts,
      Fs2,
      main = "Surface Plot",
      ylab = "",
      col = viridisLite::viridis(128)
    )
    graphics::lines(time, rep(0, M), col = "black", lwd = 1)
    graphics::lines(time, rep(-d1, M), col = "blue", lwd = 1)
    graphics::lines(time, rep(-d1 - d1a, M), col = "green", lwd = 1)
    graphics::lines(time, rep(-d1 - d1a - dl, M), col = "red", lwd = 1)
    graphics::lines(time, rep(d3, M), col = "blue", lwd = 1)
    graphics::lines(time, rep(d3 + d3a, M), col = "green", lwd = 1)
    graphics::lines(time, rep(d3 + d3a + du, M), col = "red", lwd = 1)
  }
}

#' Phase Boxplot Data
#'
#' This function constructs the phase boxplot.
#'
#' @param warp_median fdawarp object from [time_warping] of aligned data using
#'   the median.
#' @param alpha quantile value (default=.05, i.e., 95%).
#' @param kp scalar for outlier cutoff (default=1).
#'
#' @return Returns a `phbox` object containing:
#'
#' - `median_x`: median warping function
#' - `Q1`: First quartile
#' - `Q3`: Second quartile
#' - `Q1a`: First quantile based on alpha
#' - `Q3a`: Second quantile based on alpha
#' - `minn`: minimum extreme function
#' - `maxx`: maximum extreme function
#' - `outlier_index`: indexes of outlier functions
#'
#' @keywords srvf alignment boxplot
#'
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A geometric
#'   approach to visualization of variability in functional data." Journal of
#'   the American Statistical Association in press: 1-34.
#'
#' @keywords internal
phbox_data <- function(warp_median, alpha = .05, kp = 1) {
  gam <- warp_median$warping_functions

  if (warp_median$rsamps) gam <- warp_median$gams

  M <- nrow(gam)
  N <- ncol(gam)
  lambda <- 0.5

  # amplitude median
  out <- SqrtMedian(gam)
  median_x <- out$gam_median
  psi_median <- out$median
  psi <- out$psi

  # compute phase distances
  time <- seq(0, 1, length.out = M)
  v <- matrix(0, M, N)
  binsize <- mean(diff(time))
  dx <- rep(0, N)
  for (i in 1:N) {
    psi[, i] <- sqrt(gradient(gam[, i], binsize))
    v[, i] <- inv_exp_map(psi_median, psi[, i])
    dx[i] <- sqrt(trapz(time, v[, i] ^ 2))
  }
  dx_ordering <- sort(dx, index.return = TRUE)$ix
  CR_50 <- dx_ordering[1:round(N / 2)]  # 50% central region
  m <- max(dx[CR_50])  # maximal phase distance with 50% central region

  # identify phase quartiles
  angle <- matrix(0, length(CR_50), length(CR_50))
  energy <- matrix(0, length(CR_50), length(CR_50))
  for (i in 1:(length(CR_50) - 1)) {
    for (j in (i + 1):length(CR_50)) {
      q1 <- v[, CR_50[i]]
      q3 <- v[, CR_50[j]]
      if (sum(q1) > 0)
        q1 <- q1 / sqrt(trapz(time, q1 * q1))
      if (sum(q3) > 0)
        q3 <- q3 / sqrt(trapz(time, q3 * q3))
      angle[i, j] <- trapz(time, q1 * q3)
      energy[i, j] <- (1 - lambda) * (dx[CR_50[i]] / m + dx[CR_50[j]] / m) -
        lambda * (angle[i, j] + 1)
    }
  }
  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1_index <- CR_50[maxloc[1, 1]]
  Q3_index <- CR_50[maxloc[1, 2]]
  Q1 <- gam[, Q1_index]
  Q3 <- gam[, Q3_index]
  Q1_psi <- sqrt(gradient(Q1, 1 / (M - 1)))
  Q3_psi <- sqrt(gradient(Q3, 1 / (M - 1)))

  # identify phase quantiles
  dx_ordering <- sort(dx, index.return = TRUE)$ix
  CR_alpha <- dx_ordering[1:round(N * (1 - alpha))]  # (1-alpha)% central region
  m <- max(dx[CR_alpha])  # maximal phase distance with (1-alpha)% central region
  angle <- matrix(0, length(CR_alpha), length(CR_alpha))
  energy <- matrix(0, length(CR_alpha), length(CR_alpha))
  for (i in 1:(length(CR_alpha) - 1)) {
    for (j in (i + 1):length(CR_alpha)) {
      q1 <- v[, CR_alpha[i]]
      q3 <- v[, CR_alpha[j]]
      if (sum(q1) > 0)
        q1 <- q1 / sqrt(trapz(time, q1 * q1))
      if (sum(q3) > 0)
        q3 <- q3 / sqrt(trapz(time, q3 * q3))
      angle[i, j] <- trapz(time, q1 * q3)
      energy[i, j] <- (1 - lambda) * (dx[CR_alpha[i]] / m + dx[CR_alpha[j]] / m)
      - lambda * (angle[i, j] + 1)
    }
  }
  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1a_index <- CR_alpha[maxloc[1, 1]]
  Q3a_index <- CR_alpha[maxloc[1, 2]]
  Q1a <- gam[, Q1a_index]
  Q3a <- gam[, Q3a_index]
  Q1a_psi <- sqrt(gradient(Q1a, 1 / (M - 1)))
  Q3a_psi <- sqrt(gradient(Q3a, 1 / (M - 1)))

  # check quartile and quantile going same direction
  tst <- trapz(time, v[, Q1a_index] * v[, Q1_index])
  if (tst < 0) {
    Q1a <- gam[, Q3a_index]
    Q3a <- gam[, Q1a_index]
  }

  # compute phase whiskers
  IQR <- dx[Q1_index] + dx[Q3_index]
  v1 <- v[, Q1_index]
  v3 <- v[, Q3_index]
  upper_v <- v3 + kp * IQR * v3 / sqrt(trapz(time, v3 * v3))
  lower_v <- v1 + kp * IQR * v1 / sqrt(trapz(time, v1 * v1))
  upper_psi <- exp_map(psi_median, upper_v)
  lower_psi <- exp_map(psi_median, lower_v)
  upper <- cumtrapz(time, upper_psi * upper_psi)
  lower <- cumtrapz(time, lower_psi * lower_psi)

  upper_dis <- sqrt(trapz(time, (upper_v) ^ 2))
  lower_dis <- sqrt(trapz(time, (lower_v) ^ 2))
  whisker_dis <- max(c(upper_dis, lower_dis))

  # indentify phase outliers
  outlier_index <- c()
  for (i in 1:N) {
    if (dx[dx_ordering[N + 1 - i]]> whisker_dis)
      outlier_index <- c(outlier_index, dx_ordering[N + 1 - i])
    else
      break
  }

  # identify ampitude extremes
  distance_to_upper <- rep(Inf, N)
  distance_to_lower <- rep(Inf, N)
  out_50_CR <- setdiff(setdiff(1:N, CR_50), outlier_index)
  for (i in 1:length(out_50_CR)) {
    j <- out_50_CR[i]
    distance_to_upper[j] = sqrt(trapz(time, (upper_v - v[, j]) ^ 2))
    distance_to_lower[j] = sqrt(trapz(time, (lower_v - v[, j]) ^ 2))
  }
  max_index <- which.min(distance_to_upper)
  min_index <- which.min(distance_to_lower)
  min_psi <- psi[, min_index]
  max_psi <- psi[, max_index]
  minn <- gam[, min_index]
  maxx <- gam[, max_index]

  out <- list(
    median_x = median_x,
    Q1 = Q1,
    Q3 = Q3,
    Q1a = Q1a,
    Q3a = Q3a,
    minn = minn,
    maxx = maxx,
    outlier_index = outlier_index
  )
  out$time <- time
  out$psi_median <- psi_median
  out$Q1_psi <- Q1_psi
  out$Q1a_psi <- Q1a_psi
  out$Q3_psi <- Q3_psi
  out$Q3a_psi <- Q3a_psi
  out$max_psi <- max_psi
  out$min_psi <- min_psi
  out$dist <- dx
  out$Q1a_index <- Q1a_index
  out$Q1_index <- Q1_index
  out$Q3a_index <- Q3a_index
  out$Q3_index <- Q3_index

  class(out) <- 'phbox'

  out
}
