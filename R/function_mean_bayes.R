#' Bayesian Karcher Mean Calculation
#'
#' This function calculates karcher mean of functions using Bayesian method
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time sample points of functions
#' @param times factor of length of subsample points to look at (default = 5)
#' @param group (default 1:dim(f)[2])
#' @param showplot shows plots of functions (default = T)
#' @return Returns a list containing \item{distfamily}{dist matrix}
#' \item{match.matrix}{matrix of warping functions}
#' \item{position}{position}
#' \item{mu_5}{function mean}
#' \item{rtmatrix}{rtmatrix}
#' \item{sumdist}{sumdist}
#' \item{qt.fitted}{aligned srsf functions}
#' \item{estimator}{estimator}
#' \item{estimator2}{estimator2}
#' \item{regfuncs}{registered functions}
#' @keywords bayesian
#' @concept srsf alignment
#' @references Cheng, W., Dryden, I. L., and Huang, X. (2016). Bayesian
#' registration of functions and curves. Bayesian Analysis, 11(2), 447-475.
#' @examples
#' \dontrun{
#'   out <- function_mean_bayes(simu_data$f, simu_data$time)
#' }
function_mean_bayes <- function(f, time, times = 5, group = 1:dim(f)[2], showplot = TRUE){

  cut <- 5*times
  iter <- 20
  timet <- seq(0,1,length = dim(f)[1])
  m <- length(timet)-1
  n <- dim(f)[2]
  qt.matrix <- matrix(0,m,n)

  for (j in 1:n){
    qt.matrix[,j] <- Qt.matrix(f[,j],timet)
    rescale <- sqrt(m/sum((qt.matrix[,j])^2))
    qt.matrix[,j] <- rescale*qt.matrix[,j]
  }
  row <- seq(1,m,times)
  qt.fitted.matrix <- matrix(0,m,n)

  search <- NULL
  meanq <- apply(qt.matrix,1,mean)
  for (j in 1:n)
  {
    search[j] <- Enorm(qt.matrix[,j]-meanq)
  }
  position <- which.min(search)

  mu_5 <- qt.matrix[,position]
  mu.curve <- matrix(0,iter,m+1)
  dist.matrix <- matrix(0,iter+1,n)
  for (j in 1:n){dist.matrix[1,j] <- (Enorm(mu_5-qt.matrix[,j]))^2/m}
  rtmatrix <- matrix(0,m+1,n)
  match.matrix <- matrix(0,length(row)+1,n)

  i <- 1
  diffdist <- 1000
  while (i<iter & diffdist > 0.001){
    for (j in (1:n)){
      res <- dpcode(mu_5[row],mu_5,qt.matrix[,j],times,cut)
      match <- c(res$MatchIn2,m+1)
      idy <- approx(c(row,m+1),match,method="linear",xout=1:m)$y
      idy[idy>m] <- m
      scale <- sqrt(diff(match)*(1/times))
      scalevec <- rep(scale,each = times)
      extended_idy <- ((idy-1)[-(m+1)])*times+1
      qt_5.fitted <- scalevec*((res$q2LL)[extended_idy])
      qt.fitted.matrix[,j] <- qt_5.fitted
      dist.matrix[i+1,j] <- res$NDist
      rtmatrix[,j] <- c(idy,m+1)
      match.matrix[,j]<- match
    }
    diffdist <- abs(sum(dist.matrix[i+1,])-sum(dist.matrix[i,]))
    mu_5 <- apply(qt.fitted.matrix,1,mean)
    rescale <-  sqrt(m/sum((mu_5)^2))
    mu_5 <- rescale*mu_5
    i <- i+1
  }

  estimator2 <- mu_5
  karcher.res <- findkarcherinv(match.matrix,times,round=F)
  revscalevec <- karcher.res$revscalevec
  invidy <- (karcher.res$invidy)[-(m+1)]
  invidy[invidy>=m] <- m
  mu_5 <- revscalevec*(approx(seq(m),mu_5,xout=invidy)$y)
  rescale <-  sqrt(m/sum((mu_5)^2))
  estimator <- rescale*mu_5
  reg.curve <- matrix(0,m+1,n)
  for (j in 1:n) {reg.curve[,j] <- (spline(seq(0,m),f[,j],n=times*(m+1)-1)$y)[(rtmatrix[,j]-1)*times+1]}
  crossmean <- apply(reg.curve,1,mean)

  if (showplot)
  {
    plotl <- min(f)
    plotu <- max(f)
    plot(timet,reg.curve[,1],type="l",col=group[1],main="registered functions",ylab="",ylim=c(plotl-0.1*abs(plotl),plotu+0.1*abs(plotu)))
    for ( j in 2:n){lines(timet,reg.curve[,j],col=group[j],lty=2,lwd=1.2)}
    plot(timet,crossmean,type="l",col="red",
         main="Cross sectional mean",ylab="",ylim=c(plotl-0.1*abs(plotl),plotu+0.1*abs(plotu)))
    for ( j in 1:n){lines(timet,reg.curve[,j],col="grey")}
    lines(timet,crossmean,col="red")
  }
  sumdist <- apply(dist.matrix,1,sum)

  return(list(distfamily = dist.matrix, match.matrix = match.matrix, position = position,
              mu_5 = mu_5, rtmatrix = rtmatrix, sumdist = sumdist, qt.fitted = qt.fitted.matrix,
              estimator = estimator,estimator2 = estimator2, regfuncs = reg.curve) )
}
