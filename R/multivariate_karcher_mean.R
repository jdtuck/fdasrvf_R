#' Karcher Mean of Multivariate Functional Data
#'
#' Calculates Karcher mean or median of a collection of multivariate functional
#' data using the elastic square-root velocity (srvf) framework.
#'
#' @param beta \eqn{L \times M \times N} and it is
#'   interpreted as a sample of \eqn{N} \eqn{L}-dimensional curves observed on a
#'   grid of size \eqn{M}.
#' @param lambda A numeric value specifying the elasticity. Defaults to `0.0`.
#' @param maxit maximum number of iterations
#' @param ms string defining whether the Karcher mean (`"mean"`) or Karcher
#' median (`"median"`) is returned (default = `"mean"`)
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{betamean}{mean or median curve}
#' \item{type}{string indicating whether mean or median is returned}
#' \item{betan}{aligned curves}
#' \item{q}{array of srvfs}
#' \item{qn}{array of aligned srvfs}
#' \item{gam}{array of warping functions}
#' \item{E}{energy}
#' \item{qun}{cost function}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
#'    analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine
#'    Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' out <- multivariate_karcher_mean(beta[, , 1, 1:2], maxit = 2)
#' # note: use more functions, small for speed
multivariate_karcher_mean <- function (beta, lambda = 0.0, maxit = 20, ms = "mean")
{
  if(ms!="mean"&ms!="median"){warning("ms must be either \"mean\" or \"median\". ms has been set to \"mean\"",immediate. = T)}
  if(ms!="median"){ms = "mean"}

  tmp = dim(beta)
  n = tmp[1]
  T1 = tmp[2]
  time1 <- seq(0,1,length.out=T1)
  N = tmp[3]
  q = array(0, c(n, T1, N))
  for (ii in 1:N) {
    beta1 = beta[ , , ii]
    out = curve_to_q(beta1, FALSE)
    q[, , ii] = out$q
  }

  mu = q[, , 1]
  bmu = beta[, , 1]
  delta = 0.5
  tolv = 1e-04
  told = 5 * 0.001
  itr = 1
  sumd = rep(0, maxit + 1)
  sumd[1] = Inf
  betan = array(0, c(n, T1, N))
  qn = array(0, c(n, T1, N))
  normbar = rep(0, maxit + 1)
  if(ms == "median"){ #run for median only, saves memory if getting mean
    d_i = rep(0,N) #include vector for norm calculations
    v_d = array(0, c(n, T1, N)) #include array to hold v_i / d_i
  }


  cat("\nInitializing...\n")
  gam = matrix(0,T1,N)
  for (k in 1:N) {
    out = find_rotation_seed_unqiue(mu, q[, , k], 'O', FALSE, FALSE, lambda)
    gam[, k] = out$gambest
  }

  gam = t(gam)
  gamI = SqrtMeanInverse(t(gam))
  bmu = group_action_by_gamma_coord(bmu, gamI)
  mu = curve_to_q(bmu, FALSE)$q
  mu[is.nan(mu)] <- 0

  while (itr < maxit) {
    cat(sprintf("Iteration: %d\n", itr))
    mu = mu

    for (i in 1:N) {
      q1 = q[, , i]

      out = find_rotation_seed_unqiue(mu,q1, 'O', FALSE, FALSE, lambda)
      d = sqrt(trapz(time1, (mu-out$q2best)^2))

      qn[, , i] = out$q2best
      betan[, , i] = q_to_curve(out$q2best)

      if(ms == "median"){ #run for median only, saves computation time if getting mean
        d_i[i] = sqrt(innerprod_q2(qn[,,i], qn[,,i])) #calculate sqrt of norm of v_i
      }
      sumd[itr + 1] = sumd[itr + 1] + dist^2
    }

    if(ms == "median"){#run for median only
      sumv = rowSums(qn, dims = 2)
      sum_dinv = sum(1/d_i)
      vbar = sumv/sum_dinv
    }
    else{ #run for mean only
      sumv = rowSums(qn, dims = 2)
      vbar = sumv/N
    }

    normvbar[itr] = sqrt(innerprod_q2(vbar, vbar))
    normv = normvbar[itr]
    if ((sumd[itr]-sumd[itr+1]) < 0){
      break
    } else if ((normv > tolv) && (abs(sumd[itr + 1] - sumd[itr]) > told)) {
      mu = rowMeans(qn, dims = 2)
      betamean = q_to_curve(mu)
    }
    else {
      break
    }
    itr = itr + 1
  }

  ifelse(ms=="median",type<-"Karcher Median",type<-"Karcher Mean")
  return(list(beta = beta, mu = mu, type = type, betamean = betamean, v = v, q = q,
              E=normvbar[1:itr], cent = cent, len = len, len_q = len_q,
              qun = sumd[1:itr], mean_scale = mean_scale, mean_scale_q=mean_scale_q))
}
