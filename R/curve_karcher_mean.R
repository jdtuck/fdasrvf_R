#' Karcher Mean of Curves
#'
#' Calculates Karcher mean or median of a collection of curves using the elastic
#' square-root velocity (srvf) framework.
#'
#' @param beta Array of sizes \eqn{n \times T \times N} describing \eqn{N}
#' curves of dimension \eqn{n} evaluated on \eqn{T} points
#' @param mode Open (`"O"`) or Closed (`"C"`) curves
#' @param rotated Optimize over rotation (default = `TRUE`)
#' @param scale scale curves to unit length (default = `TRUE`)
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit maximum number of iterations
#' @param ms string defining whether the Karcher mean (`"mean"`) or Karcher
#' median (`"median"`) is returned (default = `"mean"`)
#' @param parallel A boolean specifying whether to run calculations in parallel.
#'   Defaults to `TRUE`.
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{beta}{centered data}
#' \item{betamean}{mean or median curve}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' \item{gam}{array of warping functions}
#' \item{cent}{centers of original curves}
#' \item{len}{length of curves}
#' \item{len_q}{length of srvfs}
#' \item{mean_scale}{mean length}
#' \item{mean_scale_q}{mean length srvf}
#' \item{E}{energy}
#' \item{qun}{cost function}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'    Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- curve_karcher_mean(beta[, , 1, 1:2], maxit = 2, parallel=FALSE)
#' # note: use more shapes, small for speed
curve_karcher_mean <- function (beta, mode = "O", rotated = TRUE, scale = TRUE,
                                lambda = 0.0, maxit = 20, ms = "mean",
                                parallel = TRUE)
{

  if (parallel) {
    cores <- max(parallel::detectCores() - 1, 1)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
  } else
    foreach::registerDoSEQ()

  if(ms!="mean"&ms!="median"){warning("ms must be either \"mean\" or \"median\". ms has been set to \"mean\"",immediate. = T)}
  if(ms!="median"){ms = "mean"}

  mean_scale = NA
  mean_scale_q = NA
  tmp = dim(beta)
  n = tmp[1]
  T1 = tmp[2]
  N = tmp[3]
  q = array(0, c(n, T1, N))
  len = rep(0, N)
  len_q = rep(0, N)
  cent = matrix(0, n, N)
  for (ii in 1:N) {
    beta1 = beta[ , , ii]
    centroid1 = calculatecentroid(beta1)
    cent[ , ii] = -1 * centroid1
    dim(centroid1) = c(length(centroid1), 1)
    beta1 = beta1 - repmat(centroid1, 1, T1)
    beta[ , , ii] = beta1
    out = curve_to_q(beta1)
    q[, , ii] = out$q
    len[ii] = out$len
    len_q[ii] = out$lenq
  }

  mu = q[, , 1]
  bmu = beta[, , 1]
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  v = array(0, c(n, T1, N))
  normvbar = rep(0, maxit + 1)
  if(ms == "median"){ #run for median only, saves memory if getting mean
    d_i = rep(0,N) #include vector for norm calculations
    v_d = array(0, c(n, T1, N)) #include array to hold v_i / d_i
  }


  cat("\nInitializing...\n")
  gam <- foreach::foreach(n = 1:N, .combine = cbind, .packages = "fdasrvf") %dopar% {
    gam_tmp <- find_rotation_seed_unqiue(mu, q[, , n], mode, rotated, TRUE, lambda)$gambest
  }

  gam = t(gam)
  gamI = SqrtMeanInverse(t(gam))
  bmu = group_action_by_gamma_coord(bmu, gamI)
  mu = curve_to_q(bmu)$q
  mu[is.nan(mu)] <- 0

  while (itr < maxit) {
    cat(sprintf("Iteration: %d\n", itr))
    mu = mu/sqrt(innerprod_q2(mu, mu))

    if (mode=="C"){
      basis = find_basis_normal(mu)
    }

    outfor <- foreach::foreach(n = 1:N, .combine = cbind, .packages='fdasrvf') %dopar% {
      out <- karcher_calc(q[, , n], mu, basis, rotated, mode, lambda, ms)
      v <- out$v
      v_d <- out$v_d
      dist <- out$dist

      list(v, v_d, dist)
    }

    v <- unlist(outfor[1, ])
    dim(v) <- c(n, T1, N)

    v_d <- unlist(outfor[2, ])
    dim(v_d) <- c(n, T1, N)

    dist <- unlist(outfor[3, ])
    dim(dist) <- c(N)

    sumd[itr + 1] = sumd[itr + 1] + sum(dist^2)

    if(ms == "median"){#run for median only
      sumv = rowSums(v_d, dims = 2)
      sum_dinv = sum(1/d_i)
      vbar = sumv/sum_dinv
    }
    else{ #run for mean only
      sumv = rowSums(v, dims = 2)
      vbar = sumv/N
    }

    normvbar[itr] = sqrt(innerprod_q2(vbar, vbar))
    normv = normvbar[itr]
    if ((sumd[itr]-sumd[itr+1]) < 0){
      break
    } else if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > told)) {
      mu = cos(delta * normvbar[itr]) * mu + sin(delta *
                                                   normvbar[itr]) * vbar/normvbar[itr]
      if (mode == "C") {
        mu = project_curve(mu)
      }
      x = q_to_curve(mu)
      a = -1 * calculatecentroid(x)
      dim(a) = c(length(a), 1)
      betamean = x + repmat(a, 1, T1)
    }
    else {
      break
    }
    itr = itr + 1
  }

  if (!scale){
    mean_scale = prod(len)^(1/length(len))
    mean_scale_q = prod(len_q)^(1/length(len))
    x = q_to_curve(mu, mean_scale_q)
    a = -1 * calculatecentroid(x)
    dim(a) = c(length(a), 1)
    betamean = x + repmat(a, 1, T1)
    mu = curve_to_q(betamean, scale)$q
  }

  if (parallel) parallel::stopCluster(cl)

  ifelse(ms=="median",type<-"Karcher Median",type<-"Karcher Mean")
  return(list(beta = beta, mu = mu, type = type, betamean = betamean, v = v, q = q,
              E=normvbar[1:itr], cent = cent, len = len, len_q = len_q,
              qun = sumd[1:itr], mean_scale = mean_scale, mean_scale_q=mean_scale_q))
}
