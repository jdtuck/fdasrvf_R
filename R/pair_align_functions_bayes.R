#' Align two functions
#'
#' This function aligns two functions using Bayesian SRSF framework. It will align f2
#' to f1
#'
#' @param f1 function 1
#' @param f2 function 2
#' @param time sample points of functions
#' @param iter number of iterations (default = 150000)
#' @param times number of subsample points to look at (default = 5)
#' @param tau (default ceil(times*.4))
#' @param powera (default 1)
#' @param showplot shows plots of functions (default = T)
#' @return Returns a list containing \item{f1}{function 1}
#' \item{f2_q}{registered function using best}
#' \item{gam_q}{warping function best}
#' \item{f2_a}{registered fucntion using mean}
#' \item{q2_a}{warping function mean}
#' @keywords srsf alignment, bayesian
#' @references Cheng, W., Dryden, I. L., & Huang, X. (2016). Bayesian registration of functions and curves. Bayesian Analysis, 11(2), 447–475.
#' @export
#' @examples
#' data("simu_data")
#' out = pair_align_functions_bayes(simu_data$f[,1], simu_data$f[,2], simu_data$time, iter=2) # use more iterations
pair_align_functions_bayes <- function(f1, f2, time, iter=15000, times = 5,
                                       tau = ceiling(times*.4), powera=1,
                                       showplot = TRUE){

  # Default setting shall work for many situations. If convergence issues arise then adjust proposal variance tau.
  if(times == 2) {warning("Small times may lead to convergence issues.")}
  burnin <- NULL
  kappa <- 1000
  thin <- 1
  cut <- 5*times
  alpha <- 1
  beta <- 0.001
  scale <- T

  timet <- seq(0,1,length=length(input1))
  qt1_5 <- Qt.matrix(input1,timet)
  qt2_5 <- Qt.matrix(input2,timet)
  p <- length(qt1_5)
  if (p%%times!=0) {stop(cat(sprintf("Number of points on q function = %d is not a multiple of times = %d.", p,times)))}
  L <- round(length(qt1_5)/times)
  row <- times*seq(0,L-1,1)+1
  if (scale){
    rescale <- sqrt(p/sum((qt1_5)^2))
    qt1_5 <- rescale*qt1_5
    rescale <- sqrt(p/sum((qt2_5)^2))
    qt2_5 <- rescale*qt2_5
  }
  res <- DP(qt1_5[row],qt1_5,qt2_5,times,cut)
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
  log.posterior <- res$log.posterior
  dist_collect <- res$dist_collect
  kappafamily <- res$kappafamily
  bestidy <- approx(c(row,p+1),best_match,method="linear",xout=1:p)$y
  bestidy[bestidy > p] <- p
  bestidy <- c(bestidy,p+1)
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
  Meanidy <- c(Meanidy,p+1)

  reg_q <- (spline(seq(0,p),input2,n=times*(p+1)-1)$y)[(bestidy-1)*times+1]
  reg_a <- (spline(seq(0,p),input2,n=times*(p+1)-1)$y)[(Meanidy-1)*times+1]

  if (showplot){

    input3 <- qtocurve(qt1_5,timet)
    input4 <- qtocurve(qt2_5,timet)
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

    plot (timet,input1,type="l",col="black",main="",ylab="Height",xlab="t")
    lines(timet,input2,col="blue")
    legend("topleft",c("function 1","function 2"),col=c("black","blue"),lty=c(1,1))
    title("Original functions")

    plot(timet,input1,type="l",col="black",main="",ylab="Height",xlab="t")
    lines(timet,reg_q,col="blue")
    legend("topleft",c("function 1","function 2*"),col=c("black","blue"),lty=c(1,1))
    title("Registration by DP estimate")

    plot(timet,input1,type="l",col="black",main="",ylab="Height",xlab="t")
    lines(timet,reg_a,col="blue")
    legend("topleft",c("function 1","function 2*"),col=c("black","blue"),lty=c(1,1))
    title("Registration by Bayesian estimate")

    traceplot(mcmc(kappafamily[burnin:iter]))
    title("Traceplot of kappa after burn-in")

    traceplot(mcmc(log.posterior[burnin:iter]))
    title("Traceplot of log posterior after burn-in")
  }

  return(list(f1 = input1, f2_q = reg_q, gam_q = (bestidy-1)/p,
              f2_a = reg_a, gam_a = (Meanidy-1)/p))
}
