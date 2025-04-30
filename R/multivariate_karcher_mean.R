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
#' @inheritParams calc_shape_dist
#' @param maxit An integer value specifying the maximum number of iterations.
#'   Defaults to `20L`.
#' @param ms A character string specifying whether the Karcher mean ("mean") or
#'   Karcher median ("median") is returned. Defaults to `"mean"`.
#' @param exact_medoid A boolean specifying whether to compute the exact medoid
#'   from the distance matrix or as the input curve closest to the pointwise
#'   mean. Defaults to `FALSE` for saving computational time.
#' @param ncores An integer value specifying the number of cores to use for
#'   parallel computation. Defaults to `1L`. The maximum number of available
#'   cores is determined by the **parallel** package. One core is always left
#'   out to avoid overloading the system.
#' @param verbose A boolean specifying whether to print the progress of the
#'  algorithm. Defaults to `FALSE`.
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
#' - `gamma`: A numeric array of shape \eqn{L \times M \times N} storing the warping
#' functions of the aligned, possibly optimally rotated and optimally scaled, input data.
#' - `R`: A numeric array of shape \eqn{L \times L \times N} storing the rotation
#' matrices of the aligned, possibly optimally rotated and optimally scaled, input data.
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
#'   analysis of elastic curves in euclidean spaces. IEEE Transactions on
#'   Pattern Analysis and Machine Intelligence, **33** (7), 1415-1428.
#'
#' @export
#'
#' @examples
#' out <- multivariate_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more functions, small for speed
multivariate_karcher_mean <- function(beta,
                                      mode = "O",
                                      alignment = TRUE,
                                      rotation = FALSE,
                                      scale = FALSE,
                                      lambda = 0.0,
                                      maxit = 20L,
                                      ms = c("mean", "median"),
                                      exact_medoid = FALSE,
                                      ncores = 1L,
                                      verbose = FALSE)
{
  if (mode == "C" && !scale)
    cli::cli_abort("Closed curves are currently handled only on the Hilbert
                   sphere. Please set `scale = TRUE`.")

  dims <- dim(beta)
  L <- dims[1] # Dimension of codomain
  M <- dims[2] # Size of the evaluation grid
  N <- dims[3] # Sample size
  ms <- rlang::arg_match(ms)

  # Compute number of cores to use
  navail <- max(parallel::detectCores() - 1, 1)

  if (ncores > navail) {
    cli::cli_alert_warning(
      "The number of requested cores ({ncores}) is larger than the number of
      available cores ({navail}). Using the maximum number of available cores..."
    )
    ncores <- navail
  }

  # Computes SRVFs
  srvfs <- lapply(1:N, \(n) curve_to_srvf(beta[ , , n], scale = scale))
  q <- array(dim = c(L, M, N))
  for (n in 1:N)
    q[ , , n] <- srvfs[[n]]$q

  # Find medoid
  if (exact_medoid) {
    out <- curve_dist(
      beta = beta,
      mode = mode,
      alignment = alignment,
      rotation = rotation,
      scale = scale,
      ncores = ncores
    )
    dists <- as.numeric(rowSums(as.matrix(out$Da)))
  }

  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else
    foreach::registerDoSEQ()

  if (!exact_medoid) {
    qmean <- rowMeans(q, dims = 2)
    n <- NULL
    dists <- foreach::foreach(n = 1:N, .combine = "c", .packages = 'fdasrvf') %dopar% {
      find_rotation_seed_unique(
        qmean, srvfs[[n]]$q,
        mode = mode,
        alignment = alignment,
        rotation = rotation,
        scale = scale
      )$d
    }
  }

  medoid_idx <- which.min(dists)

  # Initialize mean as the medoid
  qmean <- q[, , medoid_idx]
  delta <- 0.5
  tolv <- 1e-04
  told <- 5 * 0.001
  itr <- 1
  qn <- qt <- gam <- array(0, c(L, M, N))
  normbar <- rep(0, maxit)
  sumd <- rep(0, maxit + 1)
  sumd[1] <- Inf

  while (itr < maxit) {
    if (verbose)
      cli::cli_alert_info("Iteration {itr}/{maxit}...")

    if (mode == "O" || !scale)
      basis <- NULL
    else
      basis <- find_basis_normal(qmean)

    alignment_step <- foreach::foreach(
      n = 1:N,
      .combine = cbind,
      .packages = "fdasrvf") %dopar% {
        out <- find_rotation_seed_unique(
          q1 = qmean,
          q2 = q[ , , n],
          mode = mode,
          alignment = alignment,
          rotation = rotation,
          scale = scale,
          lambda = lambda
        )
        list(d = out$d, q2n = out$q2best, gam = out$gambest, R = out$Rbest)
      }

    d <- unlist(alignment_step[1, ])
    dim(d) <- N
    sumd[itr + 1] <- sum(d^2)

    qt <- unlist(alignment_step[2, ])
    dim(qt) <- c(L, M, N)

    gam <- unlist(alignment_step[3, ])
    dim(gam) <- c(M, N)

    R <- unlist(alignment_step[4, ])
    dim(R) <- c(L, L, N)

    out <- pointwise_karcher_mean(qt, qmean,
                                  basis = basis,
                                  scale = scale,
                                  ms = ms,
                                  delta = delta)

    normv <- sqrt(innerprod_q2(out$vbar, out$vbar))
    normbar[itr] <- normv

    if (sumd[itr] - sumd[itr + 1] < 0 ||
        normv < tolv ||
        abs(sumd[itr + 1] - sumd[itr]) < told)
      break

    qn <- qt
    qmean <- out$qmean
    v <- out$v

    itr <- itr + 1
  }

  len_q <- sapply(srvfs, \(x) x$qnorm)
  qmean_norm <- prod(len_q) ^ (1 / length(len_q))
  betamean <- q_to_curve(qmean)
  betamean <- betamean - calculatecentroid(betamean)

  # Compute beta2n
  betan <- array(dim = c(L, M, N))
  for (n in 1:N) {
    scl <- 1
    # if (scale)
    #   scl <- qmean_norm / len_q[n]
    betan[ , , n] <- q_to_curve(qn[ , , n], scale = scl)
    betan[ , , n] <- betan[ , , n] - calculatecentroid(betan[ , , n])
  }

  type <- ifelse(ms == "median", "Karcher Median", "Karcher Mean")

  out <- list(
            beta = beta,
            q = q,
            betan = betan,
            mu = qmean,
            qn = qn,
            gamma = gam,
            R = R,
            betamean = betamean,
            qmean = qmean,
            type = type,
            lambda = lambda,
            len_q = len_q,
            qmean_norm = qmean_norm,
            v = v,
            scale = scale,
            mode = mode,
            alignment = alignment,
            rotation = rotation,
            E = normbar[1:itr],
            qun = sumd[1:(itr + 1)]
          )

    class(out) <- "fdacurve"
    return(out)
}
