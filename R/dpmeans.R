#' DP-Means Clustering and Alignment
#'
#' This function clusters functions and aligns them using the elastic
#' square-root velocity function (SRVF) framework and the DP-means algorithm
#' (Kulis and Jordan, 2012), the small-variance limit of a Dirichlet process
#' mixture. A new cluster is opened whenever a function's elastic distance to
#' every existing template exceeds `lambda`, so the number of clusters is
#' learned from the data. Each pass assigns and aligns the functions,
#' normalizes the warps within each cluster and updates the templates as in
#' [kmeans_align()].
#'
#' @param f A numeric matrix of shape \eqn{M \times N} specifying a sample of
#'   \eqn{N} functions observed on a grid of size \eqn{M}.
#' @param time A numeric vector of length \eqn{M} specifying the grid on which
#'   the functions are evaluated.
#' @param lambda A numeric value specifying the cluster penalty, as an elastic
#'   distance. The objective minimized is the sum of squared distances to the
#'   templates plus `lambda^2` times the number of clusters. Defaults to `NA`,
#'   in which case it is chosen as the largest distance of a function to its
#'   nearest center in a k-means++ seeding with 4 centers.
#' @param lambda_align A numeric value specifying the alignment penalty passed
#'   to [optimum.reparam()]. Defaults to `0.0`.
#' @param alignment A Boolean specifying whether to perform alignment. Defaults
#'   to `TRUE`.
#' @param smooth_data A Boolean specifying whether to smooth data using a box
#'   filter. Defaults to `FALSE`.
#' @param sparam An integer value specifying the number of box filters applied.
#'   Defaults to `25L`.
#' @param parallel A Boolean specifying whether parallel mode (using
#'   [foreach::foreach()] and the **doParallel** package) should be activated.
#'   Defaults to `FALSE`.
#' @param omethod A string specifying which method should be used to solve the
#'   optimization problem that provides estimated warping functions. Choices are
#'   `"DP"` or `"RBFGS"`. Defaults to `"DP"`.
#' @param max_iter An integer value specifying the maximum number of iterations.
#'   Defaults to `20L`.
#' @param thresh A numeric value specifying the relative change of the
#'   objective below which convergence is assumed. Defaults to `0.01`.
#' @param use_verbose A Boolean specifying whether to display information about
#'   the calculations in the console. Defaults to `FALSE`.
#'
#' @return An object of class `fdakma` which is a list containing:
#'
#' - `f0`: the original functions;
#' - `q0`: the original SRSFs;
#' - `time`: the grid;
#' - `fn`: the aligned functions of each cluster in a list;
#' - `qn`: the aligned SRSFs of each cluster in a list;
#' - `gam`: the warping functions of each cluster in a list;
#' - `gamI`: the inverse Karcher mean warping function of each cluster in a
#' list;
#' - `labels`: the cluster memberships as an integer vector;
#' - `templates`: the centroids in the original functional space, an array of
#' shape \eqn{1 \times M \times K};
#' - `templates.q`: the centroids in SRSF space, of the same shape;
#' - `distances_to_center`: the elastic distance of each function to its
#' template at the last assignment;
#' - `lambda`: the cluster penalty used;
#' - `lambda_align`: the alignment penalty;
#' - `omethod`: the optimization method;
#' - `qun`: the DP-means objective at each iteration.
#'
#' @keywords srsf alignment clustering
#' @references Kulis, B., Jordan, M. I. (2012). "Revisiting k-means: new
#'   algorithms via Bayesian nonparametrics." Proceedings of the 29th
#'   International Conference on Machine Learning.
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'   May 2011. Registration of functional data using Fisher-Rao metric,
#'   arXiv:1103.3817v2.
#'
#' @export
#' @examples
#' \dontrun{
#'   out <- dpmeans_align(growth_vel$f, growth_vel$time)
#' }
dpmeans_align <- function(f, time,
                          lambda = NA,
                          lambda_align = 0.0,
                          alignment = TRUE,
                          smooth_data = FALSE,
                          sparam = 25L,
                          parallel = FALSE,
                          omethod = c("DP", "RBFGS"),
                          max_iter = 20L,
                          thresh = 0.01,
                          use_verbose = FALSE) {
  omethod <- rlang::arg_match(omethod)

  if (!is.matrix(f))
    stop("`f` must be a matrix of functions (M x N)")

  if (parallel) {
    cores <- max(parallel::detectCores() - 1, 1, na.rm = TRUE)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  if (smooth_data)
    f <- smooth.data(f, sparam = sparam)
  M <- nrow(f)
  N <- ncol(f)
  t01 <- seq(0, 1, length.out = M)

  q <- f_to_srvf(f, time)
  opts <- list(lambda_align = lambda_align, omethod = omethod,
               alignment = alignment)

  # start with one cluster whose template is the function closest to the
  # cross-sectional mean SRVF
  mnq <- rowMeans(q)
  dqq <- sqrt(apply((q - mnq) ^ 2, 2, function(v) trapz(time, v)))
  min_ind <- which.min(dqq)
  K <- 1
  templates <- f[, min_ind, drop = FALSE]
  templates.q <- q[, min_ind, drop = FALSE]

  if (is.na(lambda)) {
    if (use_verbose)
      cli::cli_alert_info("Initializing lambda")
    lambda <- dpmeans_kpp_init(f, q, time, 4, opts)
  }
  if (use_verbose)
    cli::cli_alert_info("lambda = {format(lambda, digits = 4)}")

  cluster.id <- rep(1L, N)
  dmin <- numeric(N)
  fa <- qa <- matrix(0, M, N)
  ga <- matrix(0, M, N)
  gamI_list <- list()
  qun <- numeric(max_iter)

  for (itr in 1:max_iter) {
    if (use_verbose)
      cli::cli_alert_info("Running iteration {itr}...")

    # align every function to every current template
    al <- lapply(1:K, function(k)
      dpmeans_align_to_template(f, q, time, templates[, k], templates.q[, k],
                                opts))

    # sequential DP-means pass: a function farther than lambda from every
    # template opens a new cluster, which later functions can join
    for (n in 1:N) {
      dn <- vapply(al, function(a) a$d[n], numeric(1))
      ik <- which.min(dn)
      dk <- dn[ik]
      if (dk > lambda) {
        K <- K + 1
        templates <- cbind(templates, f[, n])
        templates.q <- cbind(templates.q, q[, n])
        al[[K]] <- dpmeans_align_to_template(f, q, time, templates[, K],
                                             templates.q[, K], opts)
        ik <- K
        dk <- al[[K]]$d[n]
      }
      cluster.id[n] <- ik
      dmin[n] <- dk
    }

    # drop clusters that lost all their functions
    used <- sort(unique(cluster.id))
    cluster.id <- match(cluster.id, used)
    K <- length(used)
    templates <- templates[, used, drop = FALSE]
    templates.q <- templates.q[, used, drop = FALSE]
    al <- al[used]

    # DP-means objective at this assignment
    qun[itr] <- sum(dmin ^ 2) + lambda ^ 2 * K

    # keep each function aligned to its own template
    for (n in 1:N) {
      a <- al[[cluster.id[n]]]
      ga[, n] <- a$gam[, n]
      fa[, n] <- a$fn[, n]
      qa[, n] <- a$qn[, n]
    }

    # normalization: center the warps of each cluster
    gamI_list <- vector("list", K)
    for (k in 1:K) {
      id <- which(cluster.id == k)
      gamI <- SqrtMeanInverse(ga[, id, drop = FALSE])
      gamI_list[[k]] <- gamI
      for (n in id) {
        fa[, n] <- warp_f_gamma(fa[, n], time, gamI)
        qa[, n] <- f_to_srvf(fa[, n], time)
        ga[, n] <- stats::approx(t01, ga[, n], xout = gamI)$y
      }
    }

    # template identification
    for (k in 1:K) {
      id <- which(cluster.id == k)
      templates.q[, k] <- rowMeans(qa[, id, drop = FALSE])
      templates[, k] <- srvf_to_f(templates.q[, k], time, mean(fa[1, id]))
    }

    if (itr > 1 && abs(qun[itr] - qun[itr - 1]) < thresh * qun[itr - 1])
      break
  }

  ftmp <- qtmp <- gamtmp <- vector("list", K)
  for (k in 1:K) {
    id <- which(cluster.id == k)
    ftmp[[k]] <- fa[, id, drop = FALSE]
    qtmp[[k]] <- qa[, id, drop = FALSE]
    gamtmp[[k]] <- ga[, id, drop = FALSE]
  }

  out <- list(
    f0 = f,
    q0 = q,
    time = time,
    fn = ftmp,
    qn = qtmp,
    gam = gamtmp,
    gamI = gamI_list,
    labels = cluster.id,
    templates = array(templates, dim = c(1, M, K)),
    templates.q = array(templates.q, dim = c(1, M, K)),
    distances_to_center = dmin,
    lambda = lambda,
    lambda_align = lambda_align,
    omethod = omethod,
    qun = qun[1:itr]
  )
  class(out) <- "fdakma"
  out
}

# align every function to one template; gam, fn and qn are M x N and d is the
# elastic distance of each function to the template
dpmeans_align_to_template <- function(f, q, time, tf, tq, opts) {
  M <- nrow(f)
  N <- ncol(f)
  n <- NULL # foreach iterator; silences R CMD check global-variable note
  res <- foreach::foreach(n = 1:N, .packages = "fdasrvf") %dopar% {
    if (opts$alignment) {
      g <- optimum.reparam(
        Q1 = tq, T1 = time, Q2 = q[, n], T2 = time,
        lambda = opts$lambda_align, method = opts$omethod,
        f1o = tf[1], f2o = f[1, n]
      )
    } else {
      g <- seq(0, 1, length.out = M)
    }
    fw <- warp_f_gamma(f[, n], time, g)
    qw <- f_to_srvf(fw, time)
    list(g = g, fw = fw, qw = qw, d = sqrt(trapz(time, (qw - tq) ^ 2)))
  }
  list(
    gam = vapply(res, `[[`, numeric(M), "g"),
    fn = vapply(res, `[[`, numeric(M), "fw"),
    qn = vapply(res, `[[`, numeric(M), "qw"),
    d = vapply(res, `[[`, numeric(1), "d")
  )
}

# choose the DP-means penalty from a k-means++ seeding: pick K centers by
# sampling proportional to the squared elastic distance to the nearest chosen
# center and return the largest distance from any function to its nearest
# center (Kulis and Jordan, 2012)
dpmeans_kpp_init <- function(f, q, time, K, opts) {
  N <- ncol(q)
  D <- matrix(0, K, N)
  idx <- sample.int(N, 1)
  for (i in 1:K) {
    D[i, ] <- dpmeans_align_to_template(f, q, time, f[, idx], q[, idx],
                                        opts)$d
    if (i < K) {
      w <- apply(D[1:i, , drop = FALSE], 2, min) ^ 2
      if (sum(w) == 0) {
        D <- D[1:i, , drop = FALSE]
        break
      }
      idx <- which(stats::runif(1) * sum(w) < cumsum(w))[1]
    }
  }
  max(apply(D, 2, min))
}
