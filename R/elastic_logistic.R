elastic.logistic <- function(f, y, time, B=NULL, df=20, max_itr=20, 
                               smooth_data = FALSE, sparam = 25, parallel = FALSE,
                               cores=2){

  if (parallel){
    library(doParallel)
    cl = makeCluster(cores)
    registerDoParallel(cl)
  } else
  {
    registerDoSEQ()
  }
  
  binsize = mean(diff(time))
  eps = .Machine$double.eps
  M = nrow(f)
  N = ncol(f)
  f0 = f
  
  if (smooth_data){
    f = smooth.data(f,sparam)
  }
  
  # Create B-Spline Basis if none provided
  if (is.null(B)){
    B = bs(time, df=df, degree=4, intercept=TRUE)
  }
  Nb = ncol(B)
  
  # Compute q-function of the functional data
  tmp = gradient.spline(f,binsize,smooth_data)
  f = tmp$f
  q = tmp$g/sqrt(abs(tmp$g)+eps)
  
  gam = kronecker(matrix(1,1,N),seq(0,1,length.out=M))
  
  itr = 1
  LL = rep(0, max_itr)
  while(itr <= max_itr) {
    cat(sprintf("Iteration: r=%d\n", itr))
    # align data
    fn = matrix(0, M, N)
    qn = matrix(0, M, N)
    for (ii in 1:N){
      fn[,ii] = approx(time,f[,ii],xout=(time[length(time)]-time[1])*gam[,ii] + time[1])$y
      qn[,ii] = f_to_srvf(fn[,ii], time)
    }
    
    # Find alpha and beta using l_bfgs
    Phi = matrix(1, N, Nb+1)
    for (ii in 1:N){
      for (jj in 2:(Nb+1)){
        Phi[ii,jj] = trapz(time,qn[,ii] * B[, jj-1])
      }
    }
    b0 = rep(0,Nb+1)
    out = optim(b0, logit_loss, gr = logit_gradient, Phi, y,
                method = "L-BFGS-B", control = list(maxit=200,pgtol=1e-10))
    b = out$par
    
    alpha = b[1]
    beta = B %*% b[2:Nb+1]
    
    # compute the Loss
    LL[itr] = logit_loss(b,Phi,y)
    
    # find gamma
    gamma_new<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
      gam = logistic_warp(beta, time, q[,k], y[k])
    }
    
    if (pvecnorm(gamma-gamma_new,2) < 1e-5){
      break
    }else{
      gamma = gamma_new
    }
    
    itr = itr + 1
  }
  gamma = gamma_new
  
  if (parallel){
    stopCluster(cl)
  }
  
  return(list(alpha=alpha, beta=beta, fn=fn, qn=qn, gamma=gamma, q=q, B=B, 
              b=b[2:length(b)], Loss=LL[1:itr], mode='logistic')  )
}