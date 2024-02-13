#' K-Means Clustering and Alignment
#'
#' This function clusters functions and aligns using the elastic square-root
#' velocity function (SRVF) framework.
#'
#' @param f Either a numeric matrix or a numeric 3D array specifying the
#'   functions that need to be jointly clustered and aligned.
#'
#'   - If a matrix, it must be of shape \eqn{M \times N}. In this case, it is
#'   interpreted as a sample of \eqn{N} curves observed on a grid of size
#'   \eqn{M}.
#'   - If a 3D array, it must be of shape \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#'
#'   If this is multidimensional functional data, it is advised that
#'   `rotation==FALSE`
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   the curves are evaluated.
#' @param K An integer value specifying the number of clusters. Defaults to
#'   `1L`.
#' @param seeds An integer vector of length `K` specifying the indices of the
#'   curves in `f` which will be chosen as initial centroids. Defaults to `NULL`
#'   in which case such indices are randomly chosen.
#' @param centroid_type A string specifying the type of centroid to compute.
#'   Choices are `"mean"` or `"medoid"`. Defaults to `"mean"`.
#' @param nonempty An integer value specifying the minimum number of curves per
#'   cluster during the assignment step. Set it to a positive value to avoid the
#'   problem of empty clusters. Defaults to `0L`.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param showplot A Boolean specifying whether to show plots. Defaults to
#'   `FALSE`.
#' @param smooth_data A Boolean specifying whether to smooth data using a box
#'   filter. Defaults to `FALSE`.
#' @param sparam An integer value specifying the number of box filters applied.
#'   Defaults to `25L`.
#' @param parallel A Boolean specifying whether parallel mode (using
#'   [foreach::foreach()] and the **doParallel** package) should be activated.
#'   Defaults to `FALSE`.
#' @param alignment A Boolean specifying whether to perform alignment. Defaults
#'   to `TRUE`.
#' @param rotation A Boolean specifying whether to perform rotation. Defaults to
#'   to `FALSE`.
#' @param scale A Boolean specifying whether to scale curves to unit length. Defaults
#'   to `TRUE`.
#' @param omethod A string specifying which method should be used to solve the
#'   optimization problem that provides estimated warping functions. Choices are
#'   `"DP"` or `"RBFGS"`. Defaults to `"DP"`.
#' @param max_iter An integer value specifying the maximum number of iterations.
#'   Defaults to `50L`.
#' @param thresh A numeric value specifying a threshold on the cost function
#'   below which convergence is assumed. Defaults to `0.01`.
#' @param use_verbose A Boolean specifying whether to display information about
#'   the calculations in the console. Defaults to `FALSE`.
#'
#' @return An object of class `fdakma` which is a list containing:
#'
#' - `f0`: the original functions;
#' - `q0`: the original SRSFs;
#' - `fn`: the aligned functions as matrices or a 3D arrays of the same shape
#' than `f0` by clusters in a list;
#' - `qn`: the aligned SRSFs as matrices or a 3D arrays of the same shape
#' than `f0` separated in clusters in a list;
#' - `labels`: the cluster memberships as an integer vector;
#' - `templates`: the centroids in the original functional space;
#' - `templates.q`: the centroids in SRSF space;
#' - `distances_to_center`: A numeric vector storing the distances of each
#' observed curve to its center;
#' - `gam`: the warping functions as matrices or a 3D arrays of the same shape
#' than `f0` by clusters in a list;
#' - `qun`: cost function value.
#'
#' @keywords srsf alignment clustering
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using Fisher-Rao metric,
#'   arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative models for
#'   functional data using phase and amplitude separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references Sangalli, L. M., et al. (2010). "k-mean alignment for curve
#'   clustering." Computational Statistics & Data Analysis 54(5): 1219-1233.
#'
#' @export
#' @examples
#' \dontrun{
#'   out <- kmeans_align(growth_vel$f, growth_vel$time, K = 2)
#' }
kmeans_align <- function(f, time,
                         K = 1L,
                         seeds = NULL,
                         centroid_type = c("mean", "medoid"),
                         nonempty = 0L,
                         lambda = 0.0,
                         showplot = FALSE,
                         smooth_data = FALSE,
                         sparam = 25L,
                         parallel = FALSE,
                         alignment = TRUE,
                         rotation = FALSE,
                         scale = TRUE,
                         omethod = c("DP", "RBFGS"),
                         max_iter = 50L,
                         thresh = 0.01,
                         use_verbose = FALSE) {
  # Initialize --------------------------------------------------------------
  omethod <- rlang::arg_match(omethod)
  centroid_type <- rlang::arg_match(centroid_type)

  if (parallel) {
    cores <- max(parallel::detectCores() - 1, 1)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  # L: Codomain dimension
  # M: Grid size
  # N: Sample size
  dims <- dim(f)
  if (length(dims) == 2) {
    dim(f) <- c(1, dims)
    dims <- dim(f)
  }

  L <- dims[1]
  M <- dims[2]
  N <- dims[3]

  if (nonempty) {
    lp.ind <- c(rbind(rep(0:(N - 1), each = K), rep(N:(N + K - 1), times = N)))
    lp.beg <- seq(0, N * K * 2 - 1, 2)
    lp.rhs <- c(rep(1, N), rep(nonempty, K))
    lp.val <- rep(1, N * K * 2)
    constr.ncol <- N * K
    constr.nrow <- N + K

    constraints.coefs <- 1 * as.matrix(Matrix::sparseMatrix(
      i = lp.ind,
      p =  c(lp.beg, N * K * 2),
      dims = c(constr.nrow, constr.ncol),
      index1 = FALSE
    ))
    constraints.dense <- matrix(
      c(lp.ind + 1, rep(1:(N * K), each = 2), lp.val),
      ncol = 3
    )
    constraints.directions <- c(rep("=", N), rep(">=", K))
    constraints.rhs <- lp.rhs
  }

  if (is.null(seeds))
    template.ind <- sample(1:N, K)
  else
    template.ind <- seeds
  templates <- array(0, dim = c(L, M, K))
  for (k in 1:K)
    templates[, , k] <- f[, , template.ind[k]]
  cluster.id <- rep(0, N)
  qun <- rep(0, max_iter)

  # Convert to SRSF
  if (smooth_data) {
    for (l in 1:L)
      f[l, , ] <- smooth.data(f[l, , ], sparam = sparam)
  }

  if (L > 1){
    q <- curve_to_q(f, scale)
  } else {
    q <- f_to_srvf(f, time)
  }

  templates.q <- array(0, dim = c(L, M, K))
  for (k in 1:K)
    templates.q[, , k] <- q[, , template.ind[k]]

  # For storing final distances to center
  dtc <- rep(NA, N)

  for (itr in 1:max_iter) {
    if (use_verbose)
      cli::cli_alert_info("Running iteration {itr}...")

    if (use_verbose)
      cli::cli_alert_info("----> Alignment step")
    gam <- list()
    Dy <- matrix(0, nrow = K, ncol = N)
    qn <- list()
    fn <- list()
    fw <- matrix(nrow = L, ncol = M)
    for (k in 1:K) {
      outfor <- foreach::foreach(n = 1:N, .combine = cbind, .packages = "fdasrvf") %dopar% {
        if (alignment) {
          if (L > 1){
            out = find_rotation_seed_unique(templates.q[, , k], q[, , n],
                                            mode='O', rotation=rotation, scale=scale)
            gam_tmp = out$gambest

          } else{
            gam_tmp <- optimum.reparam(
              Q1 = templates.q[, , k], T1 = time,
              Q2 = q[, , n], T2 = time,
              lambda = lambda,
              method = omethod,
              f1o = templates[, 1, k],
              f2o = f[, 1, n]
            )
          }

        }
        else
          gam_tmp <- seq(0, 1, length.out = M)

        if (L > 1){
          fw <- group_action_by_gamma_coord(f[, , n], gam_tmp)
          qw <- curve_to_q(fw, scale)
        } else {
          fw <- warp_f_gamma(f[1, , n], time, gam_tmp)
          qw <- f_to_srvf(fw, time)
        }

        if (L > 1){
          q1dotq2 = innerprod_q2(templates.q[, , k], qw)
          if (q1dotq2 > 1){
            q1dotq2 = 1
          } else if(q1dotq2 < -1){
            q1dotq2 = -1
          }
          if (scale){
            dist = acos(q1dotq2)
          } else {
            v = templates.q[, , k]-qw
            dist = sqrt(innerprod_q2(v, v))
          }
        } else{
          dist <- trapz(time, (qw - templates.q[1, , k])^2)
          dist <- sqrt(dist)
        }

        list(gam_tmp, fw, qw, dist)
      }
      gam[[k]] <- do.call(cbind, outfor[1, ])
      f_temp <- unlist(outfor[2, ])
      dim(f_temp) <- c(L, M, N)
      q_temp <- unlist(outfor[3, ])
      dim(q_temp) <- c(L, M, N)
      qn[[k]] <- q_temp
      fn[[k]] <- f_temp
      dtil <- unlist(outfor[4, ])
      Dy[k, ] <- dtil
    }

    if (use_verbose)
      cli::cli_alert_info("----> Assignment step")

    old.cluster.id <- cluster.id
    if (!nonempty) {
      cluster.id <- apply(Dy, 2, which.min)
    } else {
      lpSolution <- lpSolve::lp(
        direction = "min",
        objective.in = c(Dy ^ 2),
        const.dir = constraints.directions,
        const.rhs = constraints.rhs,
        all.bin = TRUE,
        dense.const = constraints.dense
      )
      clusterIndicator <- matrix(lpSolution$solution, nrow = K)
      if (lpSolution$status)
        stop("lpSolve failed to find feasible solution")
      if (!all(clusterIndicator %in% c(0, 1))) {
        cat("Matrix of cluster indicators is")
        print(clusterIndicator)
        stop("Not all entries in cluster indicator matrix are binary")
      }
      cluster.id <- apply(clusterIndicator, 2, which.max)
    }

    if (use_verbose)
      cli::cli_alert_info("----> Normalisation step")

    for (k in 1:K) {
      id <- which(cluster.id == k)
      N1 <- length(id)
      if (N1 == 0) next()
      ftmp <- fn[[k]][, , id, drop = FALSE]
      gamtmp <- gam[[k]][, id, drop = FALSE]
      gamI <- SqrtMeanInverse(gamtmp)

      fw <- matrix(nrow = L, ncol = M)
      outfor <- foreach::foreach(n = 1:N1, .combine = cbind, .packages = "fdasrvf") %dopar% {

        if (L > 1){
          fw <- group_action_by_gamma_coord(ftmp[, , n], gamI)
          qw <- curve_to_q(fw, scale)
        } else {
          fw <- warp_f_gamma(ftmp[1, , n], time, gamI)
          qw <- f_to_srvf(fw, time)
        }

        gamt1 <- stats::approx(
          x = time,
          y = gamtmp[, n],
          xout = (time[M] - time[1]) * gamI + time[1]
        )$y
        list(gamt1, fw, qw)
      }
      gam[[k]][, id] <- do.call(cbind, outfor[1, ])
      f_temp <- unlist(outfor[2, ])
      dim(f_temp) <- c(L, M, N1)
      q_temp <- unlist(outfor[3, ])
      dim(q_temp) <- c(L, M, N1)
      qn[[k]][, , id] <- q_temp
      fn[[k]][, , id] <- f_temp
    }

    if (use_verbose)
      cli::cli_alert_info("----> Template identification step")

    qun.t <- rep(0, K)
    old.templates.q <- templates.q
    for (k in 1:K) {
      id <- which(cluster.id == k)

      if (length(id) == 0) {
        qun.t[k] <- NA
        next()
      }

      if (L > 1){
        if (scale == FALSE){
          if (centroid_type == "mean") {
            out <- multivariate_karcher_mean(fn[[k]][l, , id])
            templates[1, , k] <- out$betamean
            templates.q[1, , k] <- curve_to_q(out$betamean, scale)
          } else if (centroid_type == "medoid") {
            out <- multivariate_karcher_mean(fn[[k]][l, , id], ms='median')
            templates.q[1, , k] <- out$betamean
            templates[1, , k] <- curve_to_q(out$betamean, scale)
          }
        } else {
          if (centroid_type == "mean") {
            out = curve_karcher_mean(fn[[k]][l, , id], mode = "O", rotated = rotation)
            templates[1, , k] <- out$betamean
            templates.q[1, , k] <- curve_to_q(out$betamean, scale)
          } else if (centroid_type == "medoid") {
            out = curve_karcher_mean(fn[[k]][l, , id], mode = "O", rotated = rotation,
                                     ms='median')
            templates[1, , k] <- out$betamean
            templates.q[1, , k] <- curve_to_q(out$betamean, scale)
          }
        }
      } else {
        if (centroid_type == "mean") {
          templates.q[1, , k] <- rowMeans(qn[[k]][1, , id])
          templates[1, , k] <- rowMeans(fn[[k]][1, , id])
        } else if (centroid_type == "medoid") {
          idx <- which.min(Dy[k, id])
          templates.q[1, , k] <- qn[[k]][1, , id][, idx]
          templates[1, , k] <- fn[[k]][1, , id][, idx]
        }
      }

      qun.t[k] <- pvecnorm(templates.q[, , k] - old.templates.q[, , k], 2) /
        pvecnorm(old.templates.q[, , k], 2)
    }
    qun[itr] <- mean(qun.t, na.rm = TRUE)

    for (n in 1:N)
      dtc[n] <- Dy[cluster.id[n], n]

    if (qun[itr] < thresh || sum(table(cluster.id, old.cluster.id) > 0) == K)
      break
  }

  if (use_verbose)
    cli::cli_alert_info("Consolidating output...")

  ftmp <- qtmp <- gamtmp <- list()
  for (k in 1:K) {
    id <- which(cluster.id == k)
    if (length(id) == 0)
      next()
    ftmp[[k]] <- fn[[k]][, , id]
    qtmp[[k]] <- qn[[k]][, , id]
    gamtmp[[k]] <- gam[[k]][, id]
  }

  out <- list(
    f0 = f[, , ],
    q0 = q[, , ],
    time = time,
    fn = ftmp,
    qn = qtmp,
    gam = gamtmp,
    labels = cluster.id,
    templates = templates,
    templates.q = templates.q,
    distances_to_center = dtc,
    lambda = lambda,
    omethod = omethod,
    qun = qun[1:itr]
  )

  class(out) <- "fdakma"

  if (showplot) {
    idx <- unique(out$labels)
    tbl <- table(out$labels)
    cols <- lapply(idx, function(.idx) as.integer(.idx, tbl[.idx]))
    cols <- as.integer(unlist(cols))

    oldpar <- graphics::par(mfrow = c(1, L))

    if (L == 1) {
      graphics::matplot(time, out$f0, type = "l")
    } else {
      for (l in 1:L)
        graphics::matplot(time, out$f0[l, , ], type = "l", ylab = paste("Component", l))
    }
    graphics::mtext("Curves in original functional space", side = 3, line = -2, outer = TRUE)

    if (L == 1) {
      curves <- matrix(nrow = M, ncol = 0)
      for (k in 1:K)
        curves <- cbind(curves, out$fn[[k]])
      graphics::matplot(time, curves, type = "l", col = cols)
    } else {
      for (l in 1:L) {
        curves <- matrix(nrow = M, ncol = 0)
        for (k in 1:K)
          curves <- cbind(curves, out$fn[[k]][l, , ])
        graphics::matplot(time, curves, type = "l", ylab = paste("Component", l), col = cols)
      }
    }
    graphics::mtext("Aligned curves in original functional space", side = 3, line = -2, outer = TRUE)

    if (L == 1) {
      graphics::matplot(time, out$q0, type = "l")
    } else {
      for (l in 1:L)
        graphics::matplot(time, out$q0[l, , ], type = "l", ylab = paste("Component", l))
    }
    graphics::mtext("Curves in SRSF space", side = 3, line = -2, outer = TRUE)

    if (L == 1) {
      curves <- matrix(nrow = M, ncol = 0)
      for (k in 1:K)
        curves <- cbind(curves, out$qn[[k]])
      graphics::matplot(time, curves, type = "l", col = cols)
    } else {
      for (l in 1:L) {
        curves <- matrix(nrow = M, ncol = 0)
        for (k in 1:K)
          curves <- cbind(curves, out$qn[[k]][l, , ])
        graphics::matplot(time, curves, type = "l", ylab = paste("Component", l), col = cols)
      }
    }
    graphics::mtext("Aligned curves in SRSF space", side = 3, line = -2, outer = TRUE)

    graphics::par(oldpar)

    curves <- matrix(nrow = M, ncol = 0)
    for (k in 1:K)
      curves <- cbind(curves, out$gam[[k]])
    graphics::matplot(time, curves, type = "l", ylab = "warped time", col = cols)
    graphics::mtext("Estimated warping functions", side = 3, line = -2, outer = TRUE)
  }

  if (parallel) parallel::stopCluster(cl)

  out
}
