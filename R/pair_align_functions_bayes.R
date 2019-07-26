#' Align two functions
#'
#' This function aligns two functions using Bayesian SRSF framework. It will align f2
#' to f1
#'
#' @param f1 function 1
#' @param f2 function 2
#' @param timet sample points of functions
#' @param iter number of iterations (default = 15000)
#' @param times factor of length of subsample points to look at (default = 5)
#' @param tau standard deviation of Normal prior for increment (default ceil(times*.4))
#' @param powera Dirchelet prior parameter (default 1)
#' @param showplot shows plots of functions (default = T)
#' @param extrainfo T/F whether additional information is returned
#' @return Returns a list containing \item{f1}{function 1}
#' \item{f2_q}{registered function using quotient space}
#' \item{gam_q}{warping function quotient space}
#' \item{f2_a}{registered function using ambient space}
#' \item{q2_a}{warping function ambient space}
#' \item{match_collect}{posterior samples from warping function (returned if extrainfo=TRUE)}
#' \item{dist_collect}{posterior samples from the distances (returned if extrainfo=TRUE)}
#' \item{kappa_collect}{posterior samples from kappa (returned if extrainfo=TRUE)}
#' \item{log_collect}{log-likelihood of each sample (returned if extrainfo=TRUE)}
#' \item{pct_accept}{vector of acceptance ratios for the warping function (returned if extrainfo=TRUE)}
#' @keywords bayesian
#' @concept srsf alignment
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian
#' registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @export
#' @examples
#' data("simu_data")
#' out = pair_align_functions_bayes(simu_data$f[,1], simu_data$f[,2], simu_data$time)
pair_align_functions_bayes <- function(f1, f2, timet, iter=15000, times = 5,
                                       tau = ceiling(times*.4), powera=1,
                                       showplot = TRUE, extrainfo = FALSE){

  # Default setting shall work for many situations. If convergence issues arise then adjust proposal variance tau.
  if(times == 2) {warning("Small times may lead to convergence issues.")}
  burnin <- NULL
  kappa <- 1000
  thin <- 1
  cut <- 5*times
  alpha <- 1
  beta <- 0.001
  scale <- T

  qt1_5 <- f_to_srvf(f1,timet)
  qt2_5 <- f_to_srvf(f2,timet)
  p <- length(qt1_5)
  if (p%%times!=0) {
    cat(sprintf("Resampling as number of points on q function = %d is not a multiple of times = %d.\n", p,times))
    N = floor(p/times)*times
    tmp = resample.f(f1,timet,N)
    f1 = tmp$fn
    f2 = resample.f(f2,timet,N)$fn
    timet = tmp$timet
    qt1_5 <- f_to_srvf(f1,timet)
    qt2_5 <- f_to_srvf(f2,timet)
    p <- length(qt1_5)
  }
  L <- round(length(qt1_5)/times)
  row <- times*seq(0,L-1,1)+1
  if (scale){
    rescale <- sqrt(p/sum((qt1_5)^2))
    qt1_5 <- rescale*qt1_5
    rescale <- sqrt(p/sum((qt2_5)^2))
    qt2_5 <- rescale*qt2_5
  }
  res <- dpcode(qt1_5[row],qt1_5,qt2_5,times,cut)
  match <- c(res$MatchIn2,p+1)
  match_collect <- matrix(0,iter/thin,L+1)
  best_match <- match
  dist <- NULL
  dist_collect <- rep(0,iter+1)
  idy <- round(approx(c(row,p+1),match,method="linear",xout=1:p)$y)
  idy[idy > p] <- p
  scale <- sqrt(diff(match)*(1/times))
  scalevec <- rep(scale,each = times)
  dist <-(Enorm(qt1_5-scalevec*(qt2_5)[idy]))^2/p
  dist_collect[1] <- dist
  dist.min <- dist
  kappa_collect <- rep(0,iter)
  log_collect <- rep(0,iter)

  res <- simucode(iter,p,qt1_5,qt2_5, L,tau,times,kappa,alpha,beta,powera,
                  dist, dist.min, best_match, match, thin, cut)

  best_match <- res$best_match
  match_collect <- res$match_collect
  dist_min <- res$dist_min
  log_collect <- c(res$log_collect)
  dist_collect[-1] <- c(res$dist_collect)
  kappa_collect <- c(res$kappa_collect)
  bestidy <- approx(c(row,p+1),best_match,method="linear",xout=1:p)$y
  bestidy[bestidy > p] <- p
  burnin <- round(0.5*iter/thin)
  LowerP <- NULL
  UpperP <- NULL
  MeanP <- NULL
  for (i in 1:(L+1)) {
    LowerP[i] <- quantile(match_collect[burnin:(iter/thin),i],0.025)
    UpperP[i] <- quantile(match_collect[burnin:(iter/thin),i],0.975)
    MeanP[i] <- mean(match_collect[burnin:(iter/thin),i])
  }

  Meanidy <- approx(c(row,p+1),MeanP,method="linear",xout=1:p)$y
  Meanidy[Meanidy > p] <- p

  reg_q <- (spline(seq(0,p-1),f2,n=times*(p+1)-1)$y)[(bestidy-1)*times+1]
  reg_a <- (spline(seq(0,p-1),f2,n=times*(p+1)-1)$y)[(Meanidy-1)*times+1]

  if (showplot){

    input3 <- srsf_to_f(qt1_5,timet,f1[1])
    input4 <- srsf_to_f(qt2_5,timet,f2[1])
    range <- max(input3)-min(input3)
    curve1 <- input3-mean(input3)
    curve2 <- input4-mean(input4)+1*range
    plot( timet,curve1,col="black",type="l",ylim=c(min(curve1)-0.15*range,max(curve2)+0.15*range),
          main="",xlab="t", ylab="")
    lines(timet,curve2,col="blue")
    legend("topleft",c("function 1","function 2"),col=c("black","blue"),lty=c(1,1))
    SAM=1:length(best_match)
    for (n in SAM) {
      lines(c(timet[times*(n-1)+1],timet[best_match[n]]),c(curve1[times*(n-1)+1],curve2[best_match [n]]),col="red")
    }
    title("Correspondence between 2 function")

    plot ((1:(L+1)-1)/L,(best_match-1)/p,type="l",main="",xlab="t",ylab="r(t)",col="blue")
    lines((1:(L+1)-1)/L,(LowerP-1)/p, lty=2,col="red")
    lines((1:(L+1)-1)/L,(UpperP-1)/p, lty=2,col="red")
    lines((1:(L+1)-1)/L,(MeanP-1)/p,  lty=2,col="black")
    legend("topleft",c("Quotient estimate","Pointwise mean","Pointwise 95% interval"),
           col=c("blue","black","red"),lty=c(1,2,2))

    plot (timet,f1,type="l",col="black",main="",ylab="Height",xlab="t")
    lines(timet,f2,col="blue")
    legend("topleft",c("function 1","function 2"),col=c("black","blue"),lty=c(1,1))
    title("Original functions")

    plot(timet,f1,type="l",col="black",main="",ylab="Height",xlab="t")
    lines(timet,reg_q,col="blue")
    legend("topleft",c("function 1","function 2"),col=c("black","blue"),lty=c(1,1))
    title("Registration by Quotient estimate")

    plot(timet,f1,type="l",col="black",main="",ylab="Height",xlab="t")
    lines(timet,reg_a,col="blue")
    legend("topleft",c("function 1","function 2*"),col=c("black","blue"),lty=c(1,1))
    title("Registration by Bayesian estimate")

    traceplot(mcmc(kappa_collect[burnin:iter]))
    title("Traceplot of kappa after burn-in")

    traceplot(mcmc(dist_collect[(burnin:iter)+1]))
    title("Traceplot of dist after burn-in")

    traceplot(mcmc(log_collect[burnin:iter]))
    title("Traceplot of log posterior after burn-in")
  }

  retVal <- list(f1 = f1, f2_q = reg_q, gam_q = (bestidy-1)/p,
              f2_a = reg_a, gam_a = (Meanidy-1)/p)
  if (extrainfo) {
    retVal$match_collect <- match_collect
    retVal$dist_collect <- dist_collect
    retVal$kappa_collect <- kappa_collect
    retVal$log_collect <- log_collect
    retVal$pct_accept <- res$pct_accept
  }
  return(retVal)
}
