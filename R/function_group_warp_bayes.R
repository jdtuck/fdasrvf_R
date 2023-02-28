#' Bayesian Group Warping
#'
#' This function aligns a set of functions using Bayesian SRSF framework
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time sample points of functions
#' @param iter number of iterations (default = 150000)
#' @param powera Dirichlet prior parameter (default 1)
#' @param times factor of length of subsample points to look at (default = 5)
#' @param tau standard deviation of Normal prior for increment (default ceil(times*.4))
#' @param gp number of colors in plots (defaults `seq(dim(f)[2])`)
#' @param showplot shows plots of functions (default = T)
#' @return Returns a list containing \item{f0}{original functions}
#' \item{f_q}{f aligned quotient space}
#' \item{gam_q}{warping functions quotient space}
#' \item{f_a}{f aligned ambient space}
#' \item{gam_a}{warping ambient space}
#' \item{qmn}{mean srsf}
#' @keywords bayesian
#' @concept srsf alignment
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian
#'   registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @export
#' @examples
#' \dontrun{
#'   out <- function_group_warp_bayes(simu_data$f, simu_data$time)
#' }
function_group_warp_bayes <- function(f, time, iter=50000, powera=1, times=5,
                                      tau=ceiling(times*.04), gp=seq(dim(f)[2]),
                                      showplot = TRUE){

  # Default setting shall work for many situations. If convergence issues arise then adjust proposal variance tau.
  if(times==2) {warning("Small times may lead to convergence issues.")}
  burnin <- NULL
  kappa <- 1000
  alpha <- 1
  beta <-  0.001
  var.const <- 10000
  thin <- 1
  scale <- T

  m <- dim(f)[1]-1
  if (m%%times!=0) {stop(cat(sprintf("Number of points on q function = %d is not a multiple of times = %d.", m,times)))}
  timet <- seq(0,1,length=m+1)
  n <- dim(f)[2]
  qt.matrix <- matrix(0,m,n)
  qt.fitted.matrix <- matrix(0,m,n)

  for (j in 1:n){
    qt.matrix[,j] <- Qt.matrix(f[,j],timet)
    if(scale){
      rescale <- sqrt(m/sum((qt.matrix[,j])^2))
      qt.matrix[,j] <- rescale*qt.matrix[,j]
    }
  }

  row <- seq(1,m,times)
  L <- length(row)
  match.matrix <- matrix(0,L+1,n)
  best_match.matrix <- matrix(0,L+1,n)

  res.dp <- function_mean_bayes(f,times,showplot=F)
  mu_5 <- res.dp$estimator2
  match.matrix <- res.dp$match.matrix

  MAP <- mu_5
  best_match.matrix <- match.matrix
  dist.vec <- rep(100,n)
  best.vec <- dist.vec
  sumdist <- rep(0,iter)
  kappa_collect <- rep(0,iter)
  log_collect <-  rep(0,iter)
  logmax <- 0
  mu.prior <- rep(1,m)
  cov.prior <- diag(var.const,m)
  mu.q <- matrix(0,iter/thin,m)
  mu.q.standard <- matrix(0,iter/thin,m)

  burnin <- round(0.5*iter/thin)
  AVG <- length(burnin:(iter/thin))
  temp <- itercode(iter,n,m,mu_5,match.matrix,qt.matrix,qt.fitted.matrix,L,tau,
                   times,kappa,alpha,beta,powera,best.vec,dist.vec,
                   best_match.matrix,mu.prior, var.const,sumdist,thin,mu.q,
                   mu.q.standard,logmax,burnin-1, AVG)

  mu.q.standard <- temp$mu.q.standard
  mu.q <- temp$mu.q
  MAP <- temp$MAP
  best_match.matrix <- temp$best_match.matrix
  kappafamily <- temp$kappafamily
  sumdist <- temp$sumdist
  log.posterior <- temp$log.posterior
  dist <- temp$dist

  burnin <- round(0.5*iter/thin)
  mu.est <- apply(mu.q.standard[burnin:(iter/thin),],2,mean)
  rescale <- sqrt(m/sum((mu.est)^2))
  mu.est <- rescale*mu.est
  mu.est2 <- apply(mu.q[burnin:(iter/thin),],2,mean)
  rescale <- sqrt(m/sum((mu.est2)^2))
  mu.est2 <- rescale*mu.est2

  bayes_warps <- temp$bayes_warps
  gam_q <- matrix(0,m+1,n)
  gam_a <- matrix(0,m+1,n)
  f_a <- matrix(0,m+1,n)
  f_q <- matrix(0,m+1,n)

  for (t in 1:n){
    gam_q[,t] <- stats::approx(c(row,m+1),best_match.matrix[,t],method="linear",xout=1:(m+1))$y
    f_q[,t] <- (stats::spline(seq(0,m),f[,t],n=times*(m+1)-1)$y)[(gam_q[,t]-1)*times+1]
    gam_a[,t] <- stats::approx(c(row,m+1),bayes_warps[,t],method="linear", xout=1:(m+1))$y
    f_a[,t] <- (stats::spline(seq(0,m),f[,t],n=times*(m+1)-1)$y)[(gam_a[,t]-1)*times+1]
  }

  if(showplot)
  {
    coda::traceplot(coda::mcmc(log.posterior[burnin:iter]))
    graphics::title("Trace plot of log posterior after burn-in period")

    coda::traceplot(coda::mcmc(kappafamily[burnin:iter]))
    graphics::title("Trace plot of kappa after burn-in period")

    plotl <- min(f)
    plotu <- max(f)
    plot(timet,f[,1],type="l", main="",ylab="",xlab="t",
         ylim=c(plotl-0.1*abs(plotl),plotu+0.1*abs(plotu)))
    for (t in 1:n){
      graphics::lines(timet,f[,t],col=gp[t])
    }
    graphics::title("Original functions")

    plot(timet,f_q[,1],type="l",main="",ylab="",
         xlab="t",ylim=c(plotl-0.1*abs(plotl),plotu+0.1*abs(plotu)))
    for (t in 1:n){
      graphics::lines(timet,f_q[,t],col=gp[t])
    }
    graphics::title("Quotient registered functions")

    plot(timet,f_a[,1],type="l",main="",ylab="",
         xlab="t",ylim=c(plotl-0.1*abs(plotl),plotu+0.1*abs(plotu)))
    for (t in 1:n){
      graphics::lines(timet,f_a[,t],col=gp[t])
    }
    graphics::title("Bayesian registered functions")
  }

  return(list(f0 = f, f_q = f_q, gam_q = (gam_q-1)/m, f_a = f_a ,
              gam_a = (gam_a-1)/m, qmn = mu.est))
}
