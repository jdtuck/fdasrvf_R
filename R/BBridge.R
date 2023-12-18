BBridge <- function(x=0, y=0, t0=0, T=1, N=100){
  if(T<= t0) stop("wrong times")
  dt <- (T-t0)/N
  t <- seq(t0, T, length=N+1)
  X <- c(0,cumsum(stats::rnorm(N)*sqrt(dt)))
  BB <- x + X - (t-t0)/(T-t0)*(X[N+1]-y+x)
  X <- stats::ts(BB, start=t0,deltat=dt)
  return(invisible(X))
}
