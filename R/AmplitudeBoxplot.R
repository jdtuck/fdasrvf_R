#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding
#' shooting vectors and finds the median
#'
#' @param fn matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned functions with \eqn{N} samples
#' @param fmedian vector of \eqn{M} samples of the median calculated using \code{\link{time_warping}} with median
#' @param qn matrix (\eqn{N} x \eqn{M}) of \eqn{M} of aligned srvfs
#' @param qmedian vector of \eqn{M} samples of the median calculated using \code{\link{time_warping}} with median
#' @param time vector of size \eqn{N} describing the sample points
#' @param ka scalar for outlier cutoff
#' @param showplot shows plots of functions (default = T)
#' @return Returns a list containing \item{median_y}{median function}
#' \item{Q1}{First quartile}
#' \item{Q3}{Second quartile}
#' \item{minn}{minimum extreme function}
#' \item{maxx}{maximum extreme function}
#' \item{outlier_index}{indexes of outlier functions}
#' @keywords srvf alignment boxplot
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A Geometric Approach to Visualization of Variability in Functional Data." Journal of the American Statistical Association in press: 1-34.
#' @export
#' @examples
#' data("simu_warp_median")
#' data("simu_data")
#' out = AmplitudeBoxplot(simu_warp_median$fn, simu_warp_median$fmean, simu_warp_median$qn, simu_warp_median$mqn, simu_data$time, 1)
AmplitudeBoxplot <- function(fn, fmedian, qn, qmedian, time, ka, showplot=T){
  M <- nrow(fn)
  N <- ncol(fn)
  lambda <- 0.5

  # amplitude median
  median_y <- fmedian

  # compute amplitude distances
  dy = rep(0, N)
  for (i in 1:N){
    dy[i] = sqrt(trapz(time, (qmedian-qn[,i])^2))
  }
  dy_ordering <- sort(dy, index.return = T)$ix
  CR_50 <- dy_ordering[1:round(N/2)]  # 50% central region
  m <- max(dy[CR_50])  # maximal amplitude distance with 50% central region

  # identify amplitude quartiles
  angle <- matrix(0, length(CR_50), length(CR_50))
  energy <- matrix(0, length(CR_50), length(CR_50))
  for (i in 1:(length(CR_50)-1)){
    for (j in (i+1):length(CR_50)){
      q1 <- qn[,CR_50[i]] - qmedian
      q3 <- qn[,CR_50[j]] - qmedian
      q1 <- q1/sqrt(trapz(time,q1*q1))
      q3 <- q3/sqrt(trapz(time,q3*q3))
      angle[i,j] <- trapz(time,q1*q3)
      energy[i,j] <- (1-lambda) * (dy[CR_50[i]]/m + dy[CR_50[j]]/m) - lambda * (angle[i,j]+1)
    }
  }
  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1_index <- CR_50[maxloc[1,1]]
  Q3_index <- CR_50[maxloc[1,2]]
  Q1_q <- qn[,Q1_index]
  Q3_q <- qn[,Q3_index]
  Q1 <- fn[,Q1_index]
  Q3 <- fn[,Q3_index]

  # compute amplitude whiskers
  IQR <- dy[Q1_index] + dy[Q3_index]
  v1 <- Q1_q - qmedian
  v3 <- Q3_q - qmedian
  upper_q <- Q3_q + ka * IQR * v3 / sqrt(trapz(time,v3*v3))
  lower_q <- Q1_q + ka * IQR * v1 / sqrt(trapz(time,v1*v1))
  upper <- fmedian[1]+cumtrapz(time,upper_q*abs(upper_q))
  lower <- fmedian[1]+cumtrapz(time,lower_q*abs(lower_q))

  upper_dis <- sqrt(trapz(time,(upper_q-qmedian)^2))
  lower_dis <- sqrt(trapz(time,(lower_q-qmedian)^2))
  whisker_dis <- max(c(upper_dis,lower_dis))

  # indentify amplitude outliers
  outlier_index <- c()
  for (i in 1:N){
    if (dy[dy_ordering[N+1-i]]> whisker_dis){
      outlier_index <- c(outlier_index, dy_ordering[N+1-i])
    } else {
      break
    }
  }

  # identify ampitude extremes
  distance_to_upper <- rep(Inf, N)
  distance_to_lower <- rep(Inf, N)
  out_50_CR <- setdiff(setdiff(1:N, CR_50), outlier_index)
  for (i in 1:length(out_50_CR)){
    j <- out_50_CR[i]
    distance_to_upper[j] = sqrt(trapz(time,(upper_q-qn[,j])^2))
    distance_to_lower[j] = sqrt(trapz(time,(lower_q-qn[,j])^2))
  }
  max_index <- which.max(distance_to_upper)
  min_index <- which.min(distance_to_lower)
  min_q <- qn[,min_index]
  max_q <- qn[,max_index]
  minn <- fn[,min_index]
  maxx <- fn[,max_index]

  if (showplot){
    ymin <- min(c(min(fmedian),min(Q1),min(Q3),min(upper),min(lower)))
    ymax <- max(c(max(fmedian),max(Q1),max(Q3),max(upper),max(lower)))
    plot(time, fmedian, col="black",xlab="Time",main="Amplitude Boxplot", type="l", ylim=c(ymin, ymax))
    lines(time, Q1, col="blue")
    lines(time, Q3, col="green")
    lines(time, upper, col="red")
    lines(time, lower, col="magenta")

    s <- seq(0,1,length.out=100)
    Fs2 <- matrix(0,length(time), 397)
    Fs2[,1] <- (1-s[1]) * minn + s[1] * Q1
    for (j in 2:100){
      Fs2[,j] <- (1-s[j]) * minn + s[j] * Q1
      Fs2[,99+j] <- (1-s[j]) * Q1 + s[j] * fmedian
      Fs2[,198+j] <- (1-s[j]) * fmedian + s[j] * Q3
      Fs2[,297+j] <- (1-s[j]) * Q3 + s[j] * maxx
    }
    d1<-sqrt(trapz(time,(qmedian-Q1_q)^2))
    dl<-sqrt(trapz(time,(Q1_q-min_q)^2))
    d3<-sqrt(trapz(time,(qmedian-Q3_q)^2))
    du<-sqrt(trapz(time,(Q3_q-max_q)^2))
    part1<-seq(-d1-dl,-d1,length.out=100)
    part2<-seq(-d1,0,length.out=100)
    part3<-seq(0,d3,length.out=100)
    part4<-seq(d3,d3+du,length.out=100)
    allparts<-c(part1,part2[2:100],part3[2:100],part4[2:100])

    image(time, allparts, Fs2, main="Surface Plot", ylab="", col=viridis(128))

  }

  return(list(median_y=median_y,Q1=Q1,Q3=Q3,minn=minn,maxx=maxx,outlier_index=outlier_index))
}
