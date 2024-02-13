#' Karcher Mean of Multivariate Functional Data
#'
#' Calculates the Karcher mean or median of a collection of multivariate
#' functional data using the elastic square-root velocity (SRVF) framework.
#' While most of the time, the setting does not require a metric that is
#' invariant to rotation and scale, this can be achieved through the optional
#' arguments `rotation` and `scale`.
#'
#' @param beta A numeric array of shape \eqn{L \times M \times N} specifying an
#'   \eqn{N}-sample of \eqn{L}-dimensional functional data evaluated on a same
#'   grid of size \eqn{M}.
#' @param rotation A boolean specifying whether to make the metric
#'   rotation-invariant. Defaults to `FALSE`.
#' @param scale A boolean specifying whether to make the metric scale-invariant.
#'   Defaults to `FALSE`.
#' @param maxit An integer value specifying the maximum number of iterations.
#'   Defaults to `20L`.
#' @param ms A character string specifying whether the Karcher mean ("mean") or
#'   Karcher median ("median") is returned. Defaults to `"mean"`.
#' @param ncores An integer value specifying the number of cores to use for
#'   parallel computation. Defaults to `1L`. The maximum number of available
#'   cores is determined by the **parallel** package. One core is always left
#'   out to avoid overloading the system.
#'
#' @return A list with the following components:
#' - `beta`: A numeric array of shape \eqn{L \times M \times N} storing the
#' original input data.
#' - `q`: A numeric array of shape \eqn{L \times M \times N} storing the SRVFs
#' of the input data.
#' - `betan`: A numeric array of shape \eqn{L \times M \times N} storing the
#' aligned, possibly optimally rotated and optimally scaled, input data.
#' - `qn`: A numeric array of shape \eqn{L \times M \times N} storing the SRVFs
#' of the aligned, possibly optimally rotated and optimally scaled, input data.
#' - `betamean`: A numeric array of shape \eqn{L \times M} storing the Karcher
#' mean or median of the input data.
#' - `qmean`: A numeric array of shape \eqn{L \times M} storing the Karcher mean
#' or median of the SRVFs of the input data.
#' - `type`: A character string indicating whether the Karcher mean or median
#' has been returned.
#' - `E`: A numeric vector storing the energy of the Karcher mean or median at
#' each iteration.
#' - `qun`: A numeric vector storing the cost function of the Karcher mean or
#' median at each iteration.
#'
#' @keywords srvf alignment
#'
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'   analysis of elastic curves in euclidean spaces. Pattern Analysis and
#'   Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#'
#' @export
#'
#' @examples
#' out <- multivariate_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more functions, small for speed
multivariate_karcher_mean <- function(beta,
                                      rotation = FALSE,
                                      scale = FALSE,
                                      maxit = 20L,
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

  # Find medoid
  out <- curve_dist(
    beta = beta,
    mode = "O",
    rotation = rotation,
    scale = scale,
    ncores = ncores
  )
  dists <- as.numeric(rowSums(as.matrix(out$Da)))
  medoid_idx <- which.min(dists)

  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else
    foreach::registerDoSEQ()

  ms <- rlang::arg_match(ms)

  dims <- dim(beta)
  L <- dims[1] # Dimension of codomain
  M <- dims[2] # Size of the evaluation grid
  N <- dims[3] # Sample size
  q <- array(dim = c(L, M, N)) # Array for hosting SRVFs

  n <- NULL
  preprocessing_step <- foreach::foreach(n = 1:N,
                                         .combine = cbind,
                                         .packages = "fdasrvf") %dopar% {
    beta1 <- beta[ , , n]
    out <- curve_to_q(beta1, scale = scale)
    q1 <- out$q
    if (scale)
      beta1 <- beta1 / out$len
    list(q1 = q1, beta1 = beta1)
  }

  q <- unlist(preprocessing_step[1, ])
  dim(q) <- c(L, M, N)

  beta <- unlist(preprocessing_step[2, ])
  dim(beta) <- c(L, M, N)

  # Initialize mean as the medoid
  qmean <- q[, , medoid_idx]
  betamean <- beta[, , medoid_idx]
  delta <- 0.5
  tolv <- 1e-04
  told <- 5 * 0.001
  itr <- 1
  sumd <- rep(0, maxit + 1)
  sumd[1] <- Inf
  betan <- betat <- array(0, c(L, M, N))
  qn <- qt <- array(0, c(L, M, N))
  normbar <- rep(0, maxit + 1)
  if (ms == "median") {
    # run for median only, saves memory if getting mean
    d_i <- rep(0, N) #include vector for norm calculations
    v_d <- array(0, c(L, M, N)) # include array to hold v_i / d_i
  }

  while (itr < maxit) {
    cli::cli_alert_info("Iteration {itr}/{maxit}...")

    alignment_step <- foreach::foreach(
      n = 1:N,
      .combine = cbind,
      .packages = "fdasrvf") %dopar% {
        beta_i <- beta[ , , n]
        out <- calc_shape_dist(
          beta1 = betamean,
          beta2 = beta_i,
          mode = "O",
          rotation = rotation,
          scale = scale
        )
        list(d = out$d, q2n = out$q2n, beta2n = out$beta2n)
      }

    d <- unlist(alignment_step[1, ])
    dim(d) <- c(N, 1)
    sumd[itr + 1] <- sum(d^2)

    qt <- unlist(alignment_step[2, ])
    dim(qt) <- c(L, M, N)

    if (ms == "median")
      d_i[n] <- sqrt(innerprod_q2(qt[ , , n], qt[ , , n]))

    betat <- unlist(alignment_step[3, ])
    dim(betat) <- c(L, M, N)

    # AST: stopping criteria should probably be modified when scale is true
    if (ms == "median") {
      # run for median only
      sumv <- rowSums(qt, dims = 2)
      sum_dinv <- sum(1 / d_i)
      vbar <- sumv / sum_dinv
    } else {
      # run for mean only
      sumv <- rowSums(qt, dims = 2)
      vbar <- sumv / N
    }

    normbar[itr] <- sqrt(innerprod_q2(vbar, vbar))
    normv <- normbar[itr]

    if (sumd[itr] - sumd[itr + 1] < 0 ||
        normv < tolv ||
        abs(sumd[itr + 1] - sumd[itr]) < told)
      break

    qn <- qt
    betan <- betat
    qmean <- rowMeans(qn, dims = 2)
    # Should be something like this, but numerical integration is huge and it
    # requires to record value at first time point
    # betamean <- q_to_curve(qmean) + rowMeans(beta[, 1, ])
    betamean <- rowMeans(betan, dims = 2)

    itr <- itr + 1
  }

  type <- ifelse(ms == "median", "Karcher Median", "Karcher Mean")

  list(
    beta = beta,
    q = q,
    betan = betan,
    qn = qn,
    betamean = betamean,
    qmean = qmean,
    type = type,
    E = normbar[1:itr],
    qun = sumd[1:itr]
  )
}
