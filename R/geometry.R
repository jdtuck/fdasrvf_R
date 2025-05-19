exp_map <- function(psi, v, wnorm = l2_norm){
  v_norm <- wnorm(v)
  if (sum(v_norm) == 0){
    expgam <- cos(v_norm) * psi
  } else {
    expgam <- cos(v_norm) * psi + sin(v_norm) * v / v_norm
  }

  return(expgam)
}

l2_curvenorm <- function(psi, time=seq(0,1,length.out=ncol(psi))){
  sqrt(trapz(time,apply(psi^2,2,sum)))
}

#' map square root of warping function to tangent space
#'
#'
#' @param Psi vector describing psi function at center of tangent space
#' @param psi vector describing psi function to map to tangent space
#'
#' @return A numeric array of the same length as the input array `psi` storing the
#'   shooting vector of `psi`
#'
#' @keywords srvf alignment
#' @export
inv_exp_map<-function(Psi, psi){
  ip <- inner_product(Psi, psi)
  if(ip < -1){
    ip = -1
  }else if(ip > 1){
    ip = 1
  }
  theta <- acos(ip)

  if (theta < 1e-10){
    exp_inv = rep(0,length(psi))
  } else {
    exp_inv = theta / sin(theta) * (psi-cos(theta)*Psi)
  }
  return(exp_inv)
}

l2_norm<-function(psi, time=seq(0,1,length.out=length(psi))){
  l2norm <- sqrt(trapz(time,psi*psi))
  return(l2norm)
}

inner_product<-function(psi1, psi2, time=seq(0,1,length.out=length(psi1))){
  ip <- trapz(time,psi1*psi2)
  return(ip)
}

warp_q_gamma <- function(time, q, gam){
  M = length(gam)
  gam_dev = gradient(gam, 1/(M-1))
  q_tmp = stats::approx(time,q,xout=(time[length(time)]-time[1])*gam +
                          time[1])$y*sqrt(gam_dev)
  return(q_tmp)
}

randomGamma <- function(gam,num){
  out = SqrtMean(gam)
  mu = out$mu
  psi = out$psi
  vec = out$vec

  K = stats::cov(t(vec))
  out = svd(K)
  s = out$d
  U = out$u
  n = 5
  TT = nrow(vec)
  vm = rowMeans(vec)
  time <- seq(0,1,length.out=TT)

  rgam = matrix(0,num,TT)
  for (k in 1:num){
    a = stats::rnorm(n)
    v = rep(0,length(vm))
    for (i in 1:n){
      v = v + a[i]*sqrt(s[i])*U[,i]
    }
    psi <- exp_map(mu,v)

    gam0 <- cumtrapz(time,psi*psi)
    rgam[k,] = (gam0 - min(gam0))/(max(gam0)-min(gam0))
  }
  return(rgam)
}

#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding
#' shooting vectors and finds the inverse of mean
#'
#' @param gam matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
#' @return `gamI` inverse of Karcher mean warping function
#'
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' gamI <- SqrtMeanInverse(simu_warp$warping_functions)
SqrtMeanInverse <- function(gam){
  TT = nrow(gam)
  n = ncol(gam)
  eps = .Machine$double.eps
  time <- seq(0,1,length.out=TT)

  psi = matrix(0,TT,n)
  binsize <- mean(diff(time))
  for (i in 1:n){
    psi[,i] = sqrt(gradient(gam[,i],binsize))
  }

  # Find Direction
  mu = rowMeans(psi)
  stp <- .3
  maxiter = 501
  vec = matrix(0,TT,n)
  lvm = rep(0,maxiter)
  iter <- 1

  for (i in 1:n){
    vec[,i] <- inv_exp_map(mu, psi[,i])
  }
  vbar <- rowMeans(vec)
  lvm[iter] <- l2_norm(vbar)

  while (lvm[iter]>0.00000001 & iter<maxiter){
    mu <- exp_map(mu, stp*vbar)
    iter <- iter + 1
    for (i in 1:n){
      vec[,i] <- inv_exp_map(mu, psi[,i])
    }
    vbar <- rowMeans(vec)
    lvm[iter] <- l2_norm(vbar)
  }

  gam_mu = cumtrapz(time, mu*mu)
  gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))
  gamI = invertGamma(gam_mu)
  return(gamI)
}

findkarcherinv <- function(warps, times, round = F){
  m <- dim(warps)[1]
  n <- dim(warps)[2]
  psi.m <- matrix(0,m-1,n)
  for(j in 1:n){psi.m[,j]<- sqrt(diff(warps[,j])/times)}
  w <- apply(psi.m,1,mean)
  mupsi <- w/sqrt(sum(w^2/(m-1)))
  v.m <- matrix(0,m-1,n)
  check <- 1
  while(check > 0.01){
    for (i in 1:n){
      theta <- acos(sum(mupsi*psi.m[,i]/(m-1)))
      v.m[,i] <- theta/sin(theta)*(psi.m[,i]-cos(theta)*mupsi)
    }
    vbar <- apply(v.m,1,mean)
    check <- Enorm(vbar)/sqrt(m-1)
    if (check>0){
      mupsi.update <- cos(0.01*Enorm(vbar)/sqrt(m-1))*mupsi+sin(0.01*Enorm(vbar)/sqrt(m-1))*vbar/(Enorm(vbar)/sqrt(m-1))
    }
    else { mupsi.update <- cos(0.01*Enorm(vbar)/sqrt(m-1))*mupsi}
  }
  karcher.s <- 1+c(0,cumsum(mupsi.update^2)*times)
  if(round){
    invidy <- c(round(stats::approx(karcher.s,seq(1,(m-1)*times+1,times),method="linear",xout=1:((m-1)*times))$y),(m-1)*times+1)
  }
  else{
    invidy <- c((stats::approx(karcher.s,seq(1,(m-1)*times+1,times),method="linear",xout=1:((m-1)*times))$y),(m-1)*times+1)
  }
  revscalevec <- sqrt(diff(invidy))
  return(list(invidy = invidy,revscalevec = revscalevec))
}

#' map warping function to tangent space at identity
#'
#'
#' @param gam Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the warping functions
#' @param smooth Apply smoothing before gradient
#'
#' @return A numeric array of the same shape as the input array `gamma` storing the
#'   shooting vectors of `gamma` obtained via finite differences.
#'
#' @keywords srvf alignment
#' @export
gam_to_v<-function(gam, smooth=FALSE){
  if (ndims(gam) == 0){
    TT = length(gam)
    eps = .Machine$double.eps
    time <- seq(0,1,length.out=TT)
    binsize <- mean(diff(time))

    psi = rep(0,TT)
    if (smooth) {
      tmp.spline <- stats::smooth.spline(gam)
      g <- stats::predict(tmp.spline, deriv = 1)$y / binsize
      g[g<0] = 0
      psi = sqrt(g)
    } else {
        psi = sqrt(gradient(gam,binsize))
    }

    mu = rep(1,TT)
    vec <- inv_exp_map(mu, psi)

  } else {
    TT = nrow(gam)
    n = ncol(gam)
    eps = .Machine$double.eps
    time <- seq(0,1,length.out=TT)
    binsize <- mean(diff(time))

    psi = matrix(0,TT,n)
    if (smooth) {
      g <- matrix(0, TT, n)
      for (i in 1:n) {
        tmp.spline <- stats::smooth.spline(gam[,i])
        g[, i] <- stats::predict(tmp.spline, deriv = 1)$y / binsize
        g[g[,i]<0, i] = 0
        psi[,i] = sqrt(g[, i])
      }
    } else {
      for (i in 1:n){
        psi[,i] = sqrt(gradient(gam[,i],binsize))
      }
    }

    mu = rep(1,TT)
    vec = matrix(0,TT,n)
    for (i in 1:n){
      vec[,i] <- inv_exp_map(mu, psi[,i])
    }

  }

  return(vec)
}

#' map warping function to Hilbert Sphere
#'
#'
#' @param gam Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the warping functions
#' @param smooth Apply smoothing before gradient
#'
#' @return A numeric array of the same shape as the input array `gamma` storing the
#'   shooting vectors of `gamma` obtained via finite differences.
#'
#' @keywords srvf alignment
#' @export
gam_to_psi<-function(gam, smooth=FALSE){
  if (ndims(gam) == 0){
    TT = length(gam)
    eps = .Machine$double.eps
    time <- seq(0,1,length.out=TT)
    binsize <- mean(diff(time))

    psi = rep(0,TT)
    if (smooth) {
      tmp.spline <- stats::smooth.spline(gam)
      g <- stats::predict(tmp.spline, deriv = 1)$y / binsize
      g[g<0] = 0
      psi = sqrt(g)
    } else {
      psi = sqrt(gradient(gam,binsize))
    }

  } else {
    TT = nrow(gam)
    n = ncol(gam)
    eps = .Machine$double.eps
    time <- seq(0,1,length.out=TT)
    binsize <- mean(diff(time))

    psi = matrix(0,TT,n)
    if (smooth) {
      g <- matrix(0, TT, n)
      for (i in 1:n) {
        tmp.spline <- stats::smooth.spline(gam[,i])
        g[, i] <- stats::predict(tmp.spline, deriv = 1)$y / binsize
        g[g[,i]<0, i] = 0
        psi[,i] = sqrt(g[, i])
      }
    } else {
      for (i in 1:n){
        psi[,i] = sqrt(gradient(gam[,i],binsize))
      }
    }

  }

  return(psi)
}

#' map Hilbert sphere to warping function
#'
#'
#' @param psi Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the points on the Hilbert Sphere
#'
#' @return A numeric array of the same shape as the input array `psi` storing the
#'   warping function obtained via finite integration
#'
#' @keywords srvf alignment
#' @export
psi_to_gam<-function(psi){
  if (ndims(psi) == 0){
    TT = length(psi)
    time <- seq(0,1,length.out=TT)
    gam0 <- cumtrapz(time,psi*psi)
    gam <- (gam0 - min(gam0))/(max(gam0)-min(gam0))
  } else {
    TT = nrow(psi)
    n = ncol(psi)
    time <- seq(0,1,length.out=TT)

    gam = matrix(0,TT,n)
    for (i in 1:n){
      gam0 <- cumtrapz(time,psi[,i]*psi[,i])
      gam[,i] <- (gam0 - min(gam0))/(max(gam0)-min(gam0))
    }
  }
  return(gam)
}

#' map shooting vector to warping function at identity
#'
#'
#' @param v Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the shooting vectors
#'
#' @return A numeric array of the same shape as the input array `v` storing the
#'   warping functions `v`.
#'
#' @keywords srvf alignment
#' @export
v_to_gam<-function(v){
  if (ndims(v) == 0){
    TT = length(v)
    time <- seq(0,1,length.out=TT)
    mu = rep(1,TT)
    psi <- exp_map(mu,v)
    gam0 <- cumtrapz(time,psi*psi)
    gam <- (gam0 - min(gam0))/(max(gam0)-min(gam0))
  } else {
    TT = nrow(v)
    n = ncol(v)
    time <- seq(0,1,length.out=TT)

    mu = rep(1,TT)
    gam = matrix(0,TT,n)
    for (i in 1:n){
      psi <- exp_map(mu,v[,i])
      gam0 <- cumtrapz(time,psi*psi)
      gam[,i] <- (gam0 - min(gam0))/(max(gam0)-min(gam0))
    }
  }
  return(gam)
}

#' map warping function to tangent space at identity
#'
#'
#' @param gam Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the warping functions
#' @param smooth Apply smoothing before gradient
#'
#' @return A numeric array of the same shape as the input array `gamma` storing the
#'   shooting vectors of `gamma` obtained via finite differences.
#'
#' @keywords srvf alignment
#' @export
gam_to_h<-function(gam, smooth=FALSE){
  if (ndims(gam) == 0){
    TT = length(gam)
    time <- seq(0,1,length.out=TT)
    binsize <- mean(diff(time))

    psi = rep(0,TT)
    if (smooth) {
      tmp.spline <- stats::smooth.spline(gam)
      g <- stats::predict(tmp.spline, deriv = 1)$y / binsize
      g[g<0] = 0
      psi = log(g)
      h = psi - trapz(time, psi)
    } else {
        psi = log(gradient(gam,binsize))
        h = psi - trapz(time, psi)
    }

  } else {
    TT = nrow(gam)
    n = ncol(gam)
    time <- seq(0,1,length.out=TT)
    binsize <- mean(diff(time))

    h = matrix(0,TT,n)
    if (smooth) {
      g <- matrix(0, TT, n)
      for (i in 1:n) {
        tmp.spline <- stats::smooth.spline(gam[,i])
        g[, i] <- stats::predict(tmp.spline, deriv = 1)$y / binsize
        g[g[,i]<0, i] = 0
        psi = log(g[, i])
        h[, i] = psi - trapz(time, psi)
      }
    } else {
      for (i in 1:n){
        psi = log(gradient(gam[,i],binsize))
        h[, i] = psi - trapz(time, psi)
      }
    }

  }

  return(h)
}

#' map shooting vector to warping function at identity
#'
#'
#' @param h Either a numeric vector of a numeric matrix or a numeric array
#'   specifying the shooting vectors
#'
#' @return A numeric array of the same shape as the input array `h` storing the
#'   warping functions of `h`.
#'
#' @keywords srvf alignment
#' @export
h_to_gam<-function(h){
  if (ndims(h) == 0){
    TT = length(h)
    time <- seq(0,1,length.out=TT)
    gam0 <- cumtrapz(time,exp(h))
    gam0 <- gam0 / trapz(time, exp(h))
    gam <- (gam0 - min(gam0))/(max(gam0)-min(gam0))
  } else {
    TT = nrow(h)
    n = ncol(h)
    time <- seq(0,1,length.out=TT)

    gam = matrix(0,TT,n)
    for (i in 1:n){
      gam0 <- cumtrapz(time,exp(h[,i]))
      gam0 <- gam0 / trapz(time, exp(h[,i]))
      gam[,i] <- (gam0 - min(gam0))/(max(gam0)-min(gam0))
    }
  }
  return(gam)
}
