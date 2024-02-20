#' Distance Matrix Computation
#'
#' Computes the pairwise distance matrix between a set of curves using the
#' elastic shape distance as computed by [`calc_shape_dist()`].
#'
#' @param beta A numeric array of shape \eqn{L \times M \times N} specifying the
#'   set of \eqn{N} curves of length \eqn{M} in \eqn{L}-dimensional space.
#' @inheritParams calc_shape_dist
#' @param ncores An integer value specifying the number of cores to use for
#'   parallel computation. If `ncores` is greater than the number of available
#'   cores, a warning is issued and the maximum number of available cores is
#'   used. Defaults to `1L`.
#'
#' @return A list of two objects, `Da` and `Dp`, each of class `dist` containing
#'   the amplitude and phase distances, respectively.
#' @export
#'
#' @examples
#' out <- curve_dist(beta[, , 1, 1:4])
curve_dist <- function(beta,
                       mode = "O",
                       alignment = TRUE,
                       rotation = FALSE,
                       scale = FALSE,
                       include.length = FALSE,
                       ncores = 1L) {
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

  dims <- dim(beta)
  L <- dims[1]
  M <- dims[2]
  N <- dims[3]

  for (n in 1:N)
    beta[ , , n] <- curve_to_srvf(beta[, , n], scale = scale)$q

  K <- N * (N - 1) / 2

  k <- NULL
  out <- foreach::foreach(k = 0:(K - 1), .combine = cbind, .packages = "fdasrvf") %dopar% {
    # Compute indices i and j of distance matrix from linear index k
    i <- N - 2 - floor(sqrt(-8 * k + 4 * N * (N - 1) - 7) / 2.0 - 0.5)
    j <- k + i + 1 - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2

    # Increment indices as previous ones are 0-based while R expects 1-based
    beta1 <- beta[, , i + 1]
    beta2 <- beta[, , j + 1]

    out <- find_rotation_seed_unique(
      beta1, beta2,
      mode = mode,
      alignment = alignment,
      rotation = rotation,
      scale = scale
    )
    if (alignment)
      dx <- phase_distance(out$gambest)
    else
      dx <- 0
    matrix(c(out$d, dx), ncol = 1)
  }

  Da <- out[1, ]
  attributes(Da) <- NULL
  attr(Da, "Labels") <- 1:N
  attr(Da, "Size") <- N
  attr(Da, "Diag") <- FALSE
  attr(Da, "Upper") <- FALSE
  attr(Da, "call") <- match.call()
  attr(Da, "method") <- "calc_shape_dist (amplitude)"
  class(Da) <- "dist"

  Dp <- out[2, ]
  attributes(Dp) <- NULL
  attr(Dp, "Labels") <- 1:N
  attr(Dp, "Size") <- N
  attr(Dp, "Diag") <- FALSE
  attr(Dp, "Upper") <- FALSE
  attr(Dp, "call") <- match.call()
  attr(Dp, "method") <- "calc_shape_dist (phase)"
  class(Dp) <- "dist"

  list(Da = Da, Dp = Dp)
}
