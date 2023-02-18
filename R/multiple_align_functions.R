#' Group-wise function alignment to specified mean
#'
#' This function aligns a collection of functions using the elastic square-root
#' slope (srsf) framework.
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param mu vector of size \eqn{N} that f is aligned to
#' @param lambda controls the elasticity (default = 0)
#' @param pen alignment penalty (default="roughness") options are 
#' second derivative ("roughness"), geodesic distance from id ("geodesic"), and 
#' norm from id ("norm")
#' @param showplot shows plots of functions (default = T)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using [foreach()] and
#'   `doParallel` package (default=F)
#' @param omethod optimization method (DP,DP2,RBFGS,dBayes,expBayes)
#' @param MaxItr maximum number of iterations
#' @param iter bayesian number of mcmc samples (default 2000)
#' @return Returns a fdawarp object containing \item{f0}{original functions}
#' \item{fn}{aligned functions - matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples}
#' \item{qn}{aligned SRSFs - similar structure to fn}
#' \item{q0}{original SRSF - similar structure to fn}
#' \item{fmean}{function mean or median - vector of length \eqn{N}}
#' \item{mqn}{SRSF mean or median - vector of length \eqn{N}}
#' \item{gam}{warping functions - similar structure to fn}
#' \item{orig.var}{Original Variance of Functions}
#' \item{amp.var}{Amplitude Variance}
#' \item{phase.var}{Phase Variance}
#' \item{qun}{Cost Function Value}
#' @keywords srsf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
multiple_align_functions <- function(f, time, mu, lambda = 0, pen="roughness", 
                                     showplot = TRUE, smooth_data = FALSE, sparam = 25,
                                     parallel = FALSE, omethod = "DP", MaxItr = 20, iter=2000){
  if (parallel){
    cores = detectCores()-1
    cl = makeCluster(cores)
    registerDoParallel(cl)
  } else
  {
    registerDoSEQ()
  }

  cat(sprintf("lambda = %5.1f \n",lambda))

  binsize = mean(diff(time))
  eps = .Machine$double.eps
  M = nrow(f)
  N = ncol(f)
  f0 = f
  w = 0.0

  if (smooth_data){
    f = smooth.data(f,sparam)
  }

  if (showplot){
    matplot(time,f,type="l")
    title(main="Original data")
  }

  # Compute q-function of the functional data
  tmp = gradient.spline(f,binsize,smooth_data)
  f = tmp$f
  q = tmp$g/sqrt(abs(tmp$g)+eps)

  tmp = gradient.spline(mu,binsize,smooth_data)
  mf = tmp$f
  mq = tmp$g/sqrt(abs(tmp$g)+eps)
  k <- 1

  cat(sprintf("Aligning %d functions in SRSF space...\n",N))
  outfor<-foreach(k = 1:N, .combine=cbind,.packages='fdasrvf') %dopar% {
    if (omethod=="expBayes"){
      gam <- pair_align_functions_expomap(mu, c(f[,k]), time, iter=iter)$gamma
      gam <- gam$y
    } else if (omethod=="dBayes") {
      gam <- pair_align_functions_bayes(mu, f[,k], time)$gam_a
    } else {
      gam <- optimum.reparam(mq,time,q[,k],time,lambda,pen,omethod,w,mf[1],f[1,k])
    }

    gam_dev = gradient(gam,1/(M-1))
    f_temp = approx(time,f[,k],xout=(time[length(time)]-time[1])*gam +
                      time[1])$y
    q_temp = f_to_srvf(f_temp,time)
    v <- q_temp - mq
    d <- sqrt(trapz(time,v*v))
    vtil <- v/d
    dtil <- 1/d

    list(gam,gam_dev,q_temp,f_temp,vtil,dtil)
  }

  gam = unlist(outfor[1,])
  dim(gam)=c(M,N)
  gam = t(gam)
  gam_dev = unlist(outfor[2,])
  dim(gam_dev)=c(M,N)
  gam_dev = t(gam_dev)
  q_temp = unlist(outfor[3,])
  dim(q_temp)=c(M,N)
  f_temp = unlist(outfor[4,])
  dim(f_temp)=c(M,N)
  qn = q_temp
  fn = f_temp
  tmp = (1-sqrt(gam_dev))^2
  vtil = unlist(outfor[5,]);
  dim(vtil)=c(M,N)
  dtil = unlist(outfor[6,]);
  dim(dtil)=c(1,N)



  # Aligned data & stats
  q0 = q
  mean_f0 = rowMeans(f)
  std_f0 = apply(f, 1, sd)
  mean_fn = rowMeans(fn)
  std_fn = apply(fn, 1, sd)
  mqn = mq
  fmean = mean(f0[1,])+cumtrapz(time,mqn*abs(mqn))
  gam = t(gam)
  gamI = SqrtMeanInverse(gam)


  fgam = matrix(0,M,N)
  for (ii in 1:N){
    fgam[,ii] = approx(time,fmean,xout=(time[length(time)]-time[1])*gam[,ii] +
                         time[1])$y
  }
  var_fgam = apply(fgam,1,var)

  orig.var = trapz(time,std_f0^2)
  amp.var = trapz(time,std_fn^2)
  phase.var = trapz(time,var_fgam)

  out <- list(f0=f,time=time,fn=fn,qn=qn,q0=q0,fmean=fmean,mqn=mqn,gam=gam,
              orig.var=orig.var,amp.var=amp.var,phase.var=phase.var,
              qun=0,lambda=lambda,method="mean",omethod=omethod,gamI=gamI,rsamps=F)

  class(out) <- 'fdawarp'

  if (showplot){
    plot(out)
  }

  if (parallel){
    stopCluster(cl)
  }

  return(out)

}
