#' Group-wise function alignment
#'
#' This function aligns a collection of functions using the elastic square-root
#' velocity (srvf) framework.
#'
#' @param f matrix (\eqn{M} x \eqn{N}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param lambda controls the elasticity (default = 0)
#' @param showplot shows plots of functions (default = T)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using \code{\link{foreach}} and 
#'   \code{\link{doMC}} pacakge 
#' @param cores set number of cores to use with \code{\link{doMC}} (default = 8)
#' @return Returns a list containing \item{f0}{original functions}
#' \item{fn}{aligned functions}
#' \item{qn}{aligned srvfs}
#' \item{q0}{original srvfs}
#' \item{fmean}{function mean}
#' \item{mqn}{srvf mean}
#' \item{gam}{warping functions}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric, 
#'  arXiv:1103.3817v2 [math.ST].
#' @export
#' @examples
#' out = time_warping(f,time)
#' out = time_warping(f,time,showplot=F,parallel=T,cores=2)
time_warping <- function(f, time, lambda = 0, showplot = TRUE,
	smooth_data = FALSE, sparam = 25, 
	parallel = FALSE,cores=8){
	library(numDeriv)
	library(foreach)
	if (parallel){
    library(doMC)
		registerDoMC(cores=cores)
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
  
  if (smooth_data){
  	f = smooth.data(f,sparam)
  }
	
	if (showplot){
		matplot(time,f,type="l")
		title(main="Original data")
	}
	
	# Compute q-function of the functional data
	fy = gradient(f,binsize)
	q = fy/sqrt(abs(fy)+eps)
	
	cat("\nInitializing...\n")
	mnq = rowMeans(q)
	dqq = sqrt(colSums((q - matrix(mnq,ncol=N,nrow=M))^2))
	min_ind = which.min(dqq)
	mq = q[,min_ind]
	mf = f[,min_ind]

	gam<-foreach(k = 1:N, .combine=cbind) %dopar% {
	    gam_tmp = optimum.reparam(mq,time,q[,k],time)
	}
	
	gam = t(gam)
	gamI = SqrtMeanInverse(gam)
	gamI_dev = gradient(gamI, 1/(M-1))
	mf = approx(time,mf,xout=(time[length(time)]-time[1])*gamI + time[1])$y
	mq = gradient(mf,binsize)/sqrt(abs(gradient(mf,binsize))+eps)
	mq[is.nan(mq)] <- 0
	
	# Compute Mean
	cat(sprintf("Computing Karcher mean of %d functions in SRVF space...\n",N))
	MaxItr = 20
	ds = rep(0,MaxItr+2)
	ds[1] = Inf
	qun = rep(0,MaxItr)
	tmp = matrix(0,M,MaxItr+2)
	tmp[,1] = mq
	mq = tmp
	tmp = array(0,dim=c(M,N,MaxItr+2))
	tmp[,,1] = f
	f = tmp
	tmp = array(0,dim=c(M,N,MaxItr+2))
	tmp[,,1] = q
	q = tmp
	qun = rep(0,MaxItr+1)
	for (r in 1:MaxItr){
		cat(sprintf("updating step: r=%d\n", r))
		if (r == MaxItr){
			cat("maximal number of iterations is reached. \n")
		}
		
		# Matching Step
		outfor<-foreach(k = 1:N, .combine=cbind) %dopar% {
			gam = optimum.reparam(mq[,r],time,q[,k,1],time)
			gam_dev = gradient(gam,1/(M-1))
			f_temp = approx(time,f[,k,1],xout=(time[length(time)]-time[1])*gam + 
				time[1])$y
			q_temp = gradient(f_temp,binsize)/sqrt(abs(gradient(f_temp,binsize))+eps)
			list(gam,gam_dev,q_temp,f_temp)
		}
		gam = unlist(outfor[1,]);
		dim(gam)=c(M,N)
		gam = t(gam)
		gam_dev = unlist(outfor[2,]);
		dim(gam_dev)=c(M,N)
		gam_dev = t(gam_dev)
		q_temp = unlist(outfor[3,]);
		dim(q_temp)=c(M,N)
		f_temp = unlist(outfor[4,]);
		dim(f_temp)=c(M,N)
		q[,,r+1] = q_temp
    f[,,r+1] = f_temp
		tmp = (1-sqrt(gam_dev))^2
		ds_tmp  = sum(simpson(time,(matrix(mq[,r],M,N)-q[,,r+1])^2)) + 
			lambda*sum(simpson(time, t(tmp)))
    if (is.complex(ds_tmp)){
      ds[r+1] = abs(ds_tmp)
    }
  	else{
    	ds[r+1] = ds_tmp
  	}
		
		# Minimization Step
	  # compute the mean of the matched function
		mq[,r+1] = rowMeans(q[,,r+1])
		
		qun[r] = pvecnorm(mq[,r+1]-mq[,r],2)/pvecnorm(mq[,r],2)
		if (qun[r] < 1e-2 || r >=20){
			break
		}
	}
	
	if (lambda == 0){
		cat("additional run when lambda = 0\n")
		r = r+1
		outfor<-foreach(k = 1:N, .combine=cbind) %dopar% {
			gam = optimum.reparam(mq[,r],time,q[,k,1],time)
			gam_dev = gradient(gam,1/(M-1))
			list(gam,gam_dev)
		}
		gam = unlist(outfor[1,]);
		dim(gam)=c(M,N)
		gam = t(gam)
		gam_dev = unlist(outfor[2,]);
		dim(gam_dev)=c(M,N)
		gam_dev = t(gam_dev)
		
		gamI = SqrtMeanInverse(gam)
		gamI_dev = gradient(gamI, 1/(M-1))
		mq[,r+1] = approx(time,mq[,r],xout=(time[length(time)]-time[1])*gamI + 
			time[1])$y*sqrt(gamI_dev)
		
		for (k in 1:N){
			q[,k,r+1] = approx(time,q[,k,r],xout=(time[length(time)]-time[1])*gamI + 
				time[1])$y*sqrt(gamI_dev)
			f[,k,r+1] = approx(time,f[,k,r],xout=(time[length(time)]-time[1])*gamI + 
				time[1])$y
			gam[k,] = approx(time,gam[k,],xout=(time[length(time)]-time[1])*gamI + 
				time[1])$y
		}
	}
	
	# Aligned data & stats
	fn = f[,,r+1]
	qn = q[,,r+1]
	q0 = q[,,1]
	mean_f0 = rowMeans(f[,,1]);
	std_f0 = apply(f[,,1], 1, sd)
	mean_fn = rowMeans(fn)
	std_fn = apply(fn, 1, sd)
	mqn = mq[,r+1]
	fmean = mean(f0[1,])+cumtrapz(time,mqn*abs(mqn));
	gam = t(gam)
	
	if (showplot){
		matplot((0:(M-1))/(M-1),gam,type="l",main="Warping functions",xlab="Time")
		
		matplot(time,fn,type="l",main=bquote(paste("Warped Data ",lambda == 
			.(lambda))))
		
		matplot(time,cbind(mean_f0,mean_f0+std_f0,mean_f0-std_f0),type="l",lty=1,
						col=c("blue","red","green"), 
			ylab="",main=bquote(paste("Original Data: ", Mean %+-% STD)))
		legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
					 col=c('blue','red','green'),lty=1)
			
		matplot(time,cbind(mean_fn,mean_fn+std_fn,mean_fn-std_fn),type="l",lty=1,
						col=c("blue","red","green"), 
			ylab="",main=bquote(paste("Warped Data: ",lambda == .(lambda),": ",
																Mean %+-% STD)))
		legend('topright',inset=0.01,legend=c('Mean','Mean + STD', 'Mean - STD'),
					 col=c('blue','red','green'),lty=1)
		
		plot(time,fmean,type="l",col="green",main=bquote(paste(f[mean]," ", 
			lambda == .(lambda))))
	}
  
	return(list(f0=f[,,1],fn=fn,qn=qn,q0=q0,fmean=fmean,mqn=mqn,gam=gam))
  
}
