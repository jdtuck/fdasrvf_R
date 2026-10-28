#' Align New Functions to a Fitted Template
#'
#' This function aligns new functions to the template SRVF (`mqn`) of an
#' `fdawarp` object with the same `lambda`, penalty and optimization method
#' used by [time_warping()] to fit it. When the fit smoothed the data
#' (`smooth_data = TRUE`), the new functions are smoothed the same way first.
#'
#' @param x An object of class `fdawarp` typically produced by
#'   [time_warping()].
#' @param f_new A numeric matrix of shape \eqn{M \times K} of \eqn{K} new
#'   functions sampled at `x$time`, or a numeric vector of length \eqn{M} for a
#'   single function.
#'
#' @return A list with components (one column per function):
#'
#' - `fn`: \eqn{M \times K} aligned functions;
#' - `qn`: \eqn{M \times K} aligned SRVFs;
#' - `gam`: \eqn{M \times K} warping functions.
#'
#' @export
#' @examples
#' \dontrun{
#' out <- time_warping(simu_data$f, simu_data$time)
#' aligned <- align_new(out, simu_data$f[, 1:3])
#' }
align_new <- function(x, f_new) {
  if (!inherits(x, "fdawarp") || is.null(x$mqn))
    cli::cli_abort("{.arg x} must be an {.cls fdawarp} object from {.fn time_warping}.")
  align_to_template(warp_template(x), f_new)
}

# The parts of an fdawarp fit needed to align new functions to its template
warp_template <- function(x) {
  list(
    time = x$time,
    mqn = x$mqn,
    f1o = x$fmean[1],
    lambda = x$call$lambda,
    penalty_method = x$call$penalty_method,
    optim_method = x$call$optim_method,
    smooth_data = isTRUE(x$call$smooth_data),
    sparam = x$call$sparam
  )
}

align_to_template <- function(template, f_new) {
  time <- template$time
  M <- length(time)
  if (is.null(dim(f_new)))
    f_new <- matrix(f_new, ncol = 1)
  if (nrow(f_new) != M)
    cli::cli_abort("{.arg f_new} must have {M} rows (the length of {.code time}).")
  K <- ncol(f_new)
  binsize <- mean(diff(time))
  eps <- .Machine$double.eps

  # same preprocessing as time_warping()
  if (template$smooth_data)
    f_new <- smooth.data(f_new, template$sparam)
  tmp <- gradient.spline(f_new, binsize, template$smooth_data)
  f_new <- tmp$f
  q_new <- tmp$g / sqrt(abs(tmp$g) + eps)

  fn <- matrix(0, M, K)
  gam <- matrix(0, M, K)
  for (k in 1:K) {
    gam[, k] <- optimum.reparam(
      Q1 = template$mqn,
      T1 = time,
      Q2 = q_new[, k],
      T2 = time,
      lambda = template$lambda,
      pen = template$penalty_method,
      method = template$optim_method,
      f1o = template$f1o,
      f2o = f_new[1, k]
    )
    fn[, k] <- warp_f_gamma(f_new[, k], time, gam[, k])
  }
  qn <- f_to_srvf(fn, time)

  list(fn = fn, qn = qn, gam = gam)
}
