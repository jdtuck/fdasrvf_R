#elastic_regression <- function(f, y time, B=NULL, lam=0, df=20, max_itr=20, 
#                               parallel = FALSE,cores=2){
source("utility_functions.R")
source("gradient.R")
source("regression_functions.R")
source("optimum_reparam.R")
B=NULL
lam=0
df=20
max_itr=20
parallel=FALSE
cores = 8
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

# Create B-Spline Basis if none provided
if (B==NULL){
  B = bs(time, df=df, degree=4, intercept=TRUE)
}
Nb = ncol(B)

# second derivative 
Bdiff = matrix(0,M,Nb)

for (ii in 1:Nb){
  Bdiff[,ii] = gradient(gradient(B[,ii],binsize), binsize)
}

# Compute q-function of the functional data
tmp = gradient.spline(f,binsize,smooth_data)
f = tmp$f
q = tmp$g/sqrt(abs(tmp$g)+eps)

gamma = kronecker(matrix(1,1,N),seq(0,1,length.out=M))

itr = 1
SSE = rep(0, max_itr)
while(itr <= max_itr) {
  cat(sprintf("Iteration: r=%d\n", itr))
  # align data
  fn = matrix(0, M, N)
  qn = matrix(0, M, N)
  for (ii in 1:N){
    fn[,ii] = approx(time,f[,ii],xout=(time[length(time)]-time[1])*gam[,ii] + time[1])$y
    qn[,ii] = f_to_srvf(fn[,ii], time)
  }
  
  # OLS using basis
  Phi = matrix(1, N, Nb+1)
  for (ii in 1:N){
    for (jj in 2:(Nb+1)){
      Phi[ii,jj] = trapz(time,qn[,ii] * B[, jj-1])
    }
  }
  
  R = matrix(0, Nb+1, Nb+1)
  for (ii in 2:Nb+1){
    for (jj in 2:(Nb+1)){
      R[ii,jj] = trapz(time,Bdiff[,ii-1] * Bdiff[,jj-1])
    }
  }
  
  xx = t(Phi) %*% Phi
  inv_xx = solve(xx + lam * R)
  xy = t(Phi) %*% y
  b = inv_xx %*% xy
  
  alpha = b[1]
  beta = B %*% b[2:Nb+1]
  
  # compute the SSE
  int_X = rep(0, N)
  for (ii in 1:N){
    int_X[ii] = trapz(time, qn[,ii] * beta)
  }
  
  SSE[itr] = sum((y-alpha-int_X)^2)
  
  # find gamma
  gamma_new<-foreach(k = 1:N, .combine=cbind,.packages="fdasrvf") %dopar% {
    gam = regression_warp(beta, time, q[,k], y[k], alpha)
  }
  
  if (pvecnorm(gamma-gamma_new,2) < 1e-5){
    break
  }else{
    gamma = gamma_new
  }
  
  itr = itr + 1
}

# last step with centering of gam
gamI = SqrtMeanInverse(t(gam))
gamI_dev = gradient(gamI, 1/(M-1))
beta = approx(time,beta,xout=(time[length(time)]-time[1])*gamI + 
                    time[1])$y*sqrt(gamI_dev)

for (k in 1:N){
  qn[,k] = approx(time,q[,k,r],xout=(time[length(time)]-time[1])*gamI + 
                       time[1])$y*sqrt(gamI_dev)
  fn[,k] = approx(time,f[,k,r],xout=(time[length(time)]-time[1])*gamI + 
                       time[1])$y
  gam[,k] = approx(time,gam[k,],xout=(time[length(time)]-time[1])*gamI + 
                     time[1])$y
}

if (parallel){
  stopCluster(cl)
}
#  mode = 'linear'
#  return(alpha, beta, fn, qn, gamma, q, B, b[2:length(b)], SSE[1:itr], mode)  
#}