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
  K <- N * (N - 1) / 2

  out <- foreach::foreach(k = 0:(K - 1), .combine = cbind, .packages = "fdasrvf") %dopar% {
    # Compute indices i and j of distance matrix from linear index k
    i <- N - 2 - floor(sqrt(-8 * k + 4 * N * (N - 1) - 7) / 2.0 - 0.5)
    j <- k + i + 1 - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2
    # Increment indices as prevous ones are 0-based while R expects 1-based
    out <- calc_shape_dist(
      beta1 = beta[, , i + 1],
      beta2 = beta[, , j + 1],
      mode = mode,
      rotation = rotation,
      scale = scale,
      include.length = include.length
    )
    list(out$d, out$dx)
  }

  Da <- unlist(out[1, ])
  Da <- as.numeric(Da)
  attributes(Da) <- NULL
  attr(Da, "Labels") <- 1:N
  attr(Da, "Size") <- N
  attr(Da, "Diag") <- FALSE
  attr(Da, "Upper") <- FALSE
  attr(Da, "call") <- match.call()
  attr(Da, "method") <- "calc_shape_dist (amplitude)"
  class(Da) <- "dist"

  Dp <- unlist(out[2, ])
  Dp <- as.numeric(Dp)
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
