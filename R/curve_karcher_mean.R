#' Karcher Mean of Curves
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic
#' square-root velocity (SRVF) framework.
#'
#' @param beta A numeric array of shape \eqn{L \times M \times N} specifying an
#'   \eqn{N}-sample of \eqn{L}-dimensional curves evaluated on \eqn{M} points.
#' @inheritParams calc_shape_dist
#' @param rotated A boolean specifying whether to make the metric
#'   rotation-invariant. Defaults to `FALSE`.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit An integer value specifying the maximum number of iterations.
#'   Defaults to `20L`.
#' @param ms A character string specifying whether the Karcher mean ("mean") or
#'   Karcher median ("median") is returned. Defaults to `"mean"`.
#' @param ncores An integer value specifying the number of cores to use for
#'   parallel computation. Defaults to `1L`. The maximum number of available
#'   cores is determined by the **parallel** package. One core is always left
#'   out to avoid overloading the system.
#'
#' @return An object of class `fdacurve` which is a list with the following
#'   components:
#'   - beta: A numeric array of shape \eqn{L \times M \times N} specifying the
#'   \eqn{N}-sample of \eqn{L}-dimensional curves evaluated on \eqn{M} points.
#'   The curves might be slightly different from the input curves as they have
#'   been centered.
#'   - `mu`: A numeric array of shape \eqn{L \times M} specifying the Karcher
#'   mean or median of the SRVFs of the input curves.
#'   - `type`: A character string specifying whether the Karcher mean or median
#'   is returned.
#'   - `betamean`: A numeric array of shape \eqn{L \times M} specifying the
#'   Karcher mean or median of the input curves.
#'   - `v`: A numeric array of shape \eqn{L \times M \times N} specifying the
#'   shooting vectors.
#'   - `q`: A numeric array of shape \eqn{L \times M \times N} specifying the
#'   SRVFs of the input curves.
#'   - `E`: A numeric vector of shape \eqn{N} specifying XXX (TO DO)
#'   - `cent`: A numeric array of shape \eqn{L \times M} specifying the centers
#'   of the input curves.
#'   - `len`: A numeric vector of shape \eqn{N} specifying the length of the
#'   input curves.
#'   - `len_q`: A numeric vector of shape \eqn{N} specifying the length of the
#'   SRVFs of the input curves.
#'   - `qun`: A numeric vector of shape \eqn{maxit} specifying the consecutive
#'   values of the cost function.
#'   - `mean_scale`: A numeric value specifying the mean length of the input
#'   curves.
#'   - `mean_scale_q`: A numeric value specifying the mean length of the SRVFs
#'   of the input curves.
#'   - `mode`: A character string specifying the mode of the input curves.
#'   - `rotated`: A boolean specifying whether the metric is rotation-invariant.
#'   - `scale`: A boolean specifying whether the metric is scale-invariant.
#'   - `ms`: A character string specifying whether the Karcher mean or median is
#'   returned.
#'   - `lambda`: A numeric value specifying the elasticity ??? (TO DO).
#'   - `rsamps`: ??? (TO DO).
#'
#' @keywords srvf alignment
#'
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'   analysis of elastic curves in Euclidean spaces. Pattern Analysis and
#'   Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#'
#' @export
#'
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more shapes, small for speed
curve_karcher_mean <- function(beta,
                               mode = "O",
                               rotated = TRUE,
                               scale = TRUE,
                               lambda = 0.0,
                               maxit = 20,
                               ms = c("mean", "median"),
                               ncores = 1L)
{
  navail <- max(parallel::detectCores() - 1, 1)

  if (ncores > navail) {
    cli::cli_alert_warning(
      "The number of requested cores ({ncores}) is larger than the number of
      available cores ({navail}). Using the maximum number of available cores..."
    )
    ncores <- navail
  }

  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else
    foreach::registerDoSEQ()

  ms <- rlang::arg_match(ms)

  mean_scale <- NA
  mean_scale_q <- NA
  dims <- dim(beta)
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]
  q <- array(0, dim = c(L, M, N))
  len <- rep(0, N)
  len_q <- rep(0, N)
  cent <- matrix(0, nrow = L, ncol = N)

  for (n in 1:N) {
    beta1 <- beta[ , , n]
    centroid1 <- calculatecentroid(beta1)
    cent[, n] <- -1 * centroid1
    dim(centroid1) <- c(length(centroid1), 1)
    beta1 <- beta1 - repmat(centroid1, 1, M)
    beta[ , , n] <- beta1
    out <- curve_to_q(beta1, scale = TRUE)
    q[, , n] <- out$q
    len[n] <- out$len
    len_q[n] <- out$lenq
  }

  mu <- q[, , 1]
  bmu <- beta[, , 1]
  delta <- 0.5
  tolv <- 1e-04
  told <- 5 * 0.001
  itr <- 1L
  sumd <- rep(0, maxit + 1)
  sumd[1] <- Inf
  v <- array(0, dim = c(L, M, N))
  normvbar <- rep(0, maxit + 1)

  if (ms == "median") {
    # run for median only, saves memory if getting mean
    d_i <- rep(0, N) # include vector for norm calculations
    v_d <- array(0, dim = c(L, M, N)) # include array to hold v_i / d_i
  }

  cli::cli_alert_info("Initializing...")
  gam <- foreach::foreach(n = 1:N, .combine = cbind, .packages = "fdasrvf") %dopar% {
    find_rotation_seed_unique(
      q1 = mu,
      q2 = q[, , n],
      mode = mode,
      rotation = rotated,
      scale = TRUE,
      lambda = lambda
    )$gambest
  }

  gamI <- SqrtMeanInverse(gam)
  bmu <- group_action_by_gamma_coord(bmu, gamI)
  mu <- curve_to_q(bmu)$q
  mu[is.nan(mu)] <- 0

  while (itr < maxit) {
    cli::cli_alert_info("Iteration {itr}/{maxit}...")

    mu <- mu / sqrt(innerprod_q2(mu, mu))

    if (mode == "C") basis <- find_basis_normal(mu)

    outfor <- foreach::foreach(n = 1:N, .combine = cbind, .packages='fdasrvf') %dopar% {
      out <- karcher_calc(
        q1 = q[, , n],
        mu = mu,
        basis = basis,
        rotated = rotated,
        mode = mode,
        lambda = lambda,
        ms = ms
      )
      list(out$v, out$v_d, out$dist, out$d_i)
    }

    v <- unlist(outfor[1, ])
    dim(v) <- c(L, M, N)

    v_d <- unlist(outfor[2, ])
    dim(v_d) <- c(L, M, N)

    dist <- unlist(outfor[3, ])
    dim(dist) <- N

    d_i <- unlist(outfor[4, ])
    dim(d_i) <- N

    sumd[itr + 1] = sumd[itr + 1] + sum(dist^2)

    if (ms == "median") {
      # run for median only
      sumv <- rowSums(v_d, dims = 2)
      sum_dinv <- sum(1 / d_i)
      vbar <- sumv / sum_dinv
    } else {
      # run for mean only
      sumv <- rowSums(v, dims = 2)
      vbar <- sumv / N
    }

    normvbar[itr] <- sqrt(innerprod_q2(vbar, vbar))
    normv <- normvbar[itr]

    if ((sumd[itr] - sumd[itr + 1]) < 0 ||
        normv <= tolv ||
        abs(sumd[itr + 1] - sumd[itr]) < told)
      break

    mu <- cos(delta * normvbar[itr]) * mu +
      sin(delta * normvbar[itr]) * vbar / normvbar[itr]
    if (mode == "C") mu <- project_curve(mu)
    x <- q_to_curve(mu)
    a <- -1 * calculatecentroid(x)
    dim(a) <- c(length(a), 1)
    betamean <- x + repmat(a, 1, M)

    itr <- itr + 1L
  }

  if (!scale) {
    mean_scale <- prod(len) ^ (1 / length(len))
    mean_scale_q <- prod(len_q) ^ (1 / length(len))
    x <- q_to_curve(mu, mean_scale_q)
    a <- -1 * calculatecentroid(x)
    dim(a) <- c(length(a), 1)
    betamean <- x + repmat(a, 1, M)
    mu <- curve_to_q(betamean, scale)$q
  }

  type <- ifelse(ms == "median", "Karcher Median", "Karcher Mean")

  out <- list(
    beta = beta,
    mu = mu,
    type = type,
    betamean = betamean,
    v = v,
    q = q,
    E = normvbar[1:itr],
    cent = cent,
    len = len,
    len_q = len_q,
    qun = sumd[1:itr],
    mean_scale = mean_scale,
    mean_scale_q = mean_scale_q,
    mode = mode,
    rotated = rotated,
    scale = scale,
    ms = ms,
    lambda = lambda,
    rsamps = FALSE
  )

  class(out) <- "fdacurve"
  return(out)
}
