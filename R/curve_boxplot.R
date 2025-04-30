#' Curve Boxplot
#'
#' This function computes the required statistics for building up a boxplot of
#' the aligned curve data. The computed boxplot focuses on the aligned curves.
#'
#' The function returns optionally an object of class either
#' `curvebox`
#'
#' @param x An object of class `fdacurve` typically produced by [multivariate_karcher_mean()]
#' @param alpha A numeric value specifying the quantile value. Defaults to
#'   \eqn{0.05} which uses the \eqn{95\%} quantile.
#' @param range A positive numeric value specifying how far the plot whiskers
#'   extend out from the box. The whiskers extend to the most extreme data point
#'   which is no more than `range` times the interquartile range from the box.
#'   Defaults to `1.0`.
#' @param what A string specifying what the function should return. Choices are
#'   `"plot"`, `"stats"` or `"plot+stats"`. Defaults to `"plot"`.
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
#' out <- multivariate_karcher_mean(beta[, , 1, ], ms="median")
#' curve_boxplot(out, what = "stats")
#' }
curve_boxplot <- function(x,
                          alpha = 0.05,
                          range = 1.0,
                          what = c("stats", "plot","plot+stats"),
                          ...) {
  what <- rlang::arg_match(what)

  # Compute Karcher Median

  if (x$type != "Karcher Median") {
    cli::cli_alert_warning(
      "The argument {.arg x} is of class {.cls fdacurve} but has not been
      computed using the median as centroid type."
    )
    cli::cli_alert_info(
      'Rerunning {.fn multivariate_srvf_align} with {.code ms = "median"}...'
    )
    x <- multivariate_srvf_align(x$beta, mode=x$mode, rotation=x$rotation, scale=x$scale, lambda=x$lambda, ms="median")
  }

  plot_data <- curvebox_data(x, alpha = alpha, ka = range)

  if (what == "plot" || what == "plot+stats")
    boxplot(plot_data)

  if (what == "stats" || what == "plot+stats")
    return(plot_data)
}

#' @rdname boxplot.fdawarp
#' @export
boxplot.curvebox <- function(x, ...) {
  fmedian <- x$fmedian
  n = dim(fmedian)[1]
  M = dim(fmedian)[2]
  maxx <- x$maxx
  minn <- x$minn
  Q1 <- x$Q1
  Q1a <- x$Q1a
  Q3 <- x$Q3
  Q3a <- x$Q3a
  ymin <- min(c(min(fmedian), min(Q1), min(Q3), min(maxx), min(minn)))
  ymax <- max(c(max(fmedian), max(Q1), max(Q3), max(maxx), max(minn)))
  plot(fmedian[1,],
       fmedian[2,],
       col = "black",
       main = "Shape Boxplot",
       type = "l",
  )
  graphics::lines(Q1[1,], Q1[2,], col = "blue")
  graphics::lines(Q3[1,], Q3[2,], col = "blue")
  graphics::lines(Q1a[1,], Q1a[2,], col = "green")
  graphics::lines(Q3a[1,], Q3a[2,], col = "green")
  graphics::lines(maxx[1,], maxx[2,], col = "red")
  graphics::lines(minn[1,], minn[2,], col = "red")
}

#' Curve Boxplot Data
#'
#' This function constructs the amplitude boxplot.
#'
#' @param align_median object from [multivariate_karcher_mean()] of aligned curves using
#'   the median
#' @param alpha quantile value (default=.05, i.e., 95%)
#' @param ka scalar for outlier cutoff (default=1)
#'
#' @return Returns an `curvebox` object containing:
#'
#' - `median_y`: median curve
#' - `Q1`: First quartile
#' - `Q3`: Second quartile
#' - `Q1a`: First quantile based on alpha
#' - `Q3a`: Second quantile based on alpha
#' - `minn`: minimum extreme curve
#' - `maxx`: maximum extreme curve
#' - `outlier_index`: indexes of outlier curves
#' - `fmedian`: median curve
#'
#' @keywords srvf alignment boxplot
#'
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A geometric
#'   approach to visualization of variability in functional data." Journal of
#'   the American Statistical Association in press: 1-34.
#'
#' @keywords internal
curvebox_data <- function(align_median, alpha = 0.05, ka = 1) {
  fn <- align_median$betan
  median_y <- align_median$betamean
  qn <- align_median$qn
  qmedian <- align_median$mu
  scale <- align_median$scale

  if (align_median$rsamps) {
    fn <- align_median$betas
    qn <- align_median$qs
  }

  tmp = dim(fn)
  n = tmp[1]
  M = tmp[2]
  N = tmp[3]
  lambda <- 0.5

  # compute shape distances
  dy <- rep(0, N)
  for (i in 1:N){
    if (scale){
      q1dotq2 = innerprod_q2(qmedian, qn[, ,i])

      if (q1dotq2 > 1){
        q1dotq2 = 1
      } else if(q1dotq2 < -1){
        q1dotq2 = -1
      }
      dy[i] = acos(q1dotq2)
    } else {
      v = qmedian-qn[, ,i]
      dy[i] = sqrt(innerprod_q2(v, v))
    }

  }
  dy_ordering <- sort(dy, index.return = TRUE)$ix
  CR_50 <- dy_ordering[1:round(N / 2)]  # 50% central region
  m <- max(dy[CR_50])  # maximal amplitude distance with 50% central region

  # identify quartiles
  angle <- matrix(0, nrow = length(CR_50), ncol = length(CR_50))
  energy <- matrix(0, nrow = length(CR_50), ncol = length(CR_50))

  for (i in 1:(length(CR_50) - 1)) {
    for (j in (i + 1):length(CR_50)) {
      q1 <- qn[, , CR_50[i]] - qmedian
      q3 <- qn[, , CR_50[j]] - qmedian
      q1 <- q1 / sqrt(innerprod_q2(q1, q1))
      q3 <- q3 / sqrt(innerprod_q2(q3, q3))
      angle[i, j] <- innerprod_q2(q1, q3)
      energy[i, j] <- (1 - lambda) * (dy[CR_50[i]] / m + dy[CR_50[j]] / m) -
        lambda * (angle[i, j] + 1)
    }
  }

  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1_index <- CR_50[maxloc[1, 1]]
  Q3_index <- CR_50[maxloc[1, 2]]
  Q1_q <- qn[, , Q1_index]
  Q3_q <- qn[, , Q3_index]
  Q1 <- fn[, , Q1_index]
  Q3 <- fn[, , Q3_index]

  # identify amplitude quantile
  dy_ordering <- sort(dy, index.return = TRUE)$ix
  CR_alpha <- dy_ordering[1:round(N * (1 - alpha))]  # (1-alpha)% central region
  m <- max(dy[CR_alpha])  # maximal amplitude distance with (1-alpha)% central region
  angle <- matrix(0, nrow = length(CR_alpha), ncol = length(CR_alpha))
  energy <- matrix(0, nrow = length(CR_alpha), ncol = length(CR_alpha))

  for (i in 1:(length(CR_alpha) - 1)) {
    for (j in (i + 1):length(CR_alpha)) {
      q1 <- qn[, , CR_alpha[i]] - qmedian
      q3 <- qn[, , CR_alpha[j]] - qmedian
      q1 <- q1 / sqrt(innerprod_q2(q1, q1))
      q3 <- q3 / sqrt(innerprod_q2(q3, q3))
      angle[i, j] <- innerprod_q2(q1, q3)
      energy[i, j] <- (1 - lambda) * (dy[CR_alpha[i]] / m + dy[CR_alpha[j]] / m) -
        lambda * (angle[i, j] + 1)
    }
  }

  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1a_index <- CR_alpha[maxloc[1, 1]]
  Q3a_index <- CR_alpha[maxloc[1, 2]]
  Q1a_q <- qn[, , Q1a_index]
  Q3a_q <- qn[, , Q3a_index]
  Q1a <- fn[, , Q1a_index]
  Q3a <- fn[, , Q3a_index]

  # compute amplitude whiskers
  IQR <- dy[Q1_index] + dy[Q3_index]
  v1 <- Q1_q - qmedian
  v3 <- Q3_q - qmedian
  upper_q <- Q3_q + ka * IQR * v3 / sqrt(innerprod_q2(v3, v3))
  lower_q <- Q1_q + ka * IQR * v1 / sqrt(innerprod_q2(v1, v1))
  upper <- q_to_curve(upper_q)
  lower <- q_to_curve(lower_q)

  if (scale){
    q1dotq2 = innerprod_q2(upper_q, qmedian)

    if (q1dotq2 > 1){
      q1dotq2 = 1
    } else if(q1dotq2 < -1){
      q1dotq2 = -1
    }
    upper_dis = acos(q1dotq2)
  } else {
    v = upper_q-qmedian
    upper_dis = sqrt(innerprod_q2(v, v))
  }

  if (scale){
    q1dotq2 = innerprod_q2(lower_q, qmedian)

    if (q1dotq2 > 1){
      q1dotq2 = 1
    } else if(q1dotq2 < -1){
      q1dotq2 = -1
    }
    lower_dis = acos(q1dotq2)
  } else {
    v = lower_q-qmedian
    lower_dis = sqrt(innerprod_q2(v, v))
  }

  whisker_dis <- max(c(upper_dis, lower_dis))

  # identify shape outliers
  outlier_index <- c()
  for (i in 1:N) {
    if (dy[dy_ordering[N + 1 - i]] > whisker_dis)
      outlier_index <- c(outlier_index, dy_ordering[N + 1 - i])
    else
      break
  }

  # identify shape extremes
  distance_to_upper <- rep(Inf, N)
  distance_to_lower <- rep(Inf, N)
  out_50_CR <- setdiff(setdiff(1:N, CR_50), outlier_index)

  for (i in 1:length(out_50_CR)) {
    j <- out_50_CR[i]

    if (scale){
      q1dotq2 = innerprod_q2(upper_q, qn[, ,j])

      if (q1dotq2 > 1){
        q1dotq2 = 1
      } else if(q1dotq2 < -1){
        q1dotq2 = -1
      }
      distance_to_upper[j] = acos(q1dotq2)
    } else {
      v = upper_q-qn[, ,j]
      distance_to_upper[j] = sqrt(innerprod_q2(v, v))
    }

    if (scale){
      q1dotq2 = innerprod_q2(lower_q, qn[, ,j])

      if (q1dotq2 > 1){
        q1dotq2 = 1
      } else if(q1dotq2 < -1){
        q1dotq2 = -1
      }
      distance_to_lower[j] = acos(q1dotq2)
    } else {
      v = upper_q-qn[, ,j]
      distance_to_lower[j] = sqrt(innerprod_q2(v, v))
    }
  }

  max_index <- which.min(distance_to_upper)
  min_index <- which.min(distance_to_lower)
  min_q <- qn[, , min_index]
  max_q <- qn[, , max_index]
  minn <- fn[, , min_index]
  maxx <- fn[, , max_index]

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

  class(out) <- 'curvebox'

  out
}
