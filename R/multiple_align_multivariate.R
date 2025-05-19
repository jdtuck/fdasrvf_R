#' Group-wise multivariate function alignment to specified mean
#'
#' This function aligns a collection of functions using the elastic square-root
#' velocity (srvf) framework.
#'
#' @param beta A numeric array of shape \eqn{L \times M \times N} specifying an
#'   \eqn{N}-sample of \eqn{L}-dimensional functional data evaluated on a same
#'   grid of size \eqn{M}.
#' @param mu array of size \eqn{L \times M} that f is aligned to
#' @inheritParams calc_shape_dist
#' @param ncores An integer value specifying the number of cores to use for
#'   parallel computation. Defaults to `1L`. The maximum number of available
#'   cores is determined by the **parallel** package. One core is always left
#'   out to avoid overloading the system.
#' @param verbose verbose printing (default TRUE)
#' @return Returns a fdacurve object containing:
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
#' @keywords srvf alignment
#'
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'   analysis of elastic curves in euclidean spaces. IEEE Transactions on
#'   Pattern Analysis and Machine Intelligence, **33** (7), 1415-1428.
#' @export
multiple_align_multivariate <- function(beta,
                                        mu,
                                        mode = "O",
                                        alignment = TRUE,
                                        rotation = FALSE,
                                        scale = FALSE,
                                        lambda = 0.0,
                                        ncores = 1L,
                                        verbose = TRUE) {
  if (mode == "C" && !scale)
    cli::cli_abort(
      "Closed curves are currently handled only on the Hilbert
                   sphere. Please set `scale = TRUE`."
    )

  dims <- dim(beta)
  L <- dims[1] # Dimension of codomain
  M <- dims[2] # Size of the evaluation grid
  N <- dims[3] # Sample size

  # Computes SRVFs
  srvfs <- lapply(1:N, \(n) curve_to_srvf(beta[, , n], scale = scale))
  q <- array(dim = c(L, M, N))
  scales <- rep(NA, N)
  for (n in 1:N){
    scales[n] <- srvfs[[n]]$qnorm
    q[, , n] <- srvfs[[n]]$q
  }

  mq = curve_to_srvf(mu, scale = scale)$q
  k <- 1

  if (verbose) {
    cli::cli_alert_info(sprintf("Aligning %d functions in SRVF space...\n", N))
  }

  if (mode == "O" || !scale)
    basis <- NULL
  else
    basis <- find_basis_normal(mq)

  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else
    foreach::registerDoSEQ()

  alignment_step <- foreach::foreach(n = 1:N,
                                     .combine = cbind,
                                     .packages = "fdasrvf") %dopar% {
                                       out <- find_rotation_seed_unique(
                                         q1 = mq,
                                         q2 = q[, , n],
                                         mode = mode,
                                         alignment = alignment,
                                         rotation = rotation,
                                         scale = scale,
                                         lambda = lambda
                                       )
                                       list(
                                         d = out$d,
                                         q2n = out$q2best,
                                         gam = out$gambest,
                                         R = out$Rbest
                                       )
                                     }
  d <- unlist(alignment_step[1, ])
  dim(d) <- N
  sumd <- sum(d^2)

  qn <- unlist(alignment_step[2, ])
  dim(qn) <- c(L, M, N)

  gam <- unlist(alignment_step[3, ])
  dim(gam) <- c(M, N)

  R <- unlist(alignment_step[4, ])
  dim(R) <- c(L, L, N)

  len_q <- sapply(srvfs, \(x) x$qnorm)
  qmean_norm <- prod(len_q)^(1 / length(len_q))
  betamean <- mu - calculatecentroid(mu)

  # Compute beta2n
  betan <- array(dim = c(L, M, N))
  for (n in 1:N) {
    scl <- 1
    betan[, , n] <- q_to_curve(qn[, , n], scale = scl)
    betan[, , n] <- betan[, , n] - calculatecentroid(betan[, , n])
  }

  # compute shooting vectors
  v <- array(dim = c(L, M, N))
  for (n in 1:N) {
    w <- inverse_exponential_map(qn[, , n], mq, scale = scale)

    if (is.null(basis))
      v[, , n] <- w
    else
      v[, , n] <- project_tangent(w, mq, basis)
  }


  type <- "Karcher Mean"

  out <- list(
    beta = beta,
    q = q,
    betan = betan,
    mu = mq,
    qn = qn,
    gamma = gam,
    R = R,
    betamean = betamean,
    qmean = mq,
    type = type,
    lambda = lambda,
    len_q = len_q,
    qmean_norm = qmean_norm,
    v = v,
    scale = scale,
    scales = scales,
    mode = mode,
    alignment = alignment,
    rotation = rotation,
    qun = sumd
  )

  class(out) <- "fdacurve"

  return(out)
}
