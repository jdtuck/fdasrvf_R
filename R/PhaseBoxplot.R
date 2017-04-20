#' Phase Boxplot
#'
#' This function constructs the amplitude boxplot
#'
#' @param gam matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
#' @param kp scalar for outlier cutoff
#' @param showplot shows plots of functions (default = T)
#' @return Returns a list containing \item{median_x}{median warping function}
#' \item{Q1}{First quartile}
#' \item{Q3}{Second quartile}
#' \item{minn}{minimum extreme function}
#' \item{maxx}{maximum extreme function}
#' \item{outlier_index}{indexes of outlier functions}
#' @keywords srvf alignment boxplot
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A Geometric Approach to Visualization
#' of Variability in Functional Data." Journal of the American Statistical Association in press: 1-34.
#' @export
#' @examples
#' data("simu_warp_median")
#' out = PhaseBoxplot(simu_warp_median$gam, 1)
PhaseBoxplot <- function(gam, kp, showplot=T){
  M <- nrow(gam)
  N <- ncol(gam)
  lambda <- 0.5

  # amplitude median
  out <- SqrtMedian(gam)
  median_x <- out$gam_median
  psi_median <- out$median
  psi <- out$psi

  # compute phase distances
  time <- seq(0,1,length.out=M)
  v <- matrix(0,M,N)
  binsize <- mean(diff(time))
  dx = rep(0, N)
  for (i in 1:N){
    psi[,i] <- sqrt(gradient(gam[,i],binsize))
    v[,i] <- inv_exp_map(psi_median, psi[,i])
    dx[i] = sqrt(trapz(time, v[,i]^2))
  }
  dx_ordering <- sort(dx, index.return = T)$ix
  CR_50 <- dx_ordering[1:round(N/2)]  # 50% central region
  m <- max(dx[CR_50])  # maximal phase distance with 50% central region

  # identify phase quartiles
  angle <- matrix(0, length(CR_50), length(CR_50))
  energy <- matrix(0, length(CR_50), length(CR_50))
  for (i in 1:(length(CR_50)-1)){
    for (j in (i+1):length(CR_50)){
      q1 <- v[,CR_50[i]]
      q3 <- v[,CR_50[j]]
      q1 <- q1/sqrt(trapz(time,q1*q1))
      q3 <- q3/sqrt(trapz(time,q3*q3))
      angle[i,j] <- trapz(time,q1*q3)
      energy[i,j] <- (1-lambda) * (dx[CR_50[i]]/m + dx[CR_50[j]]/m) - lambda * (angle[i,j]+1)
    }
  }
  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1_index <- CR_50[maxloc[1,1]]
  Q3_index <- CR_50[maxloc[1,2]]
  Q1 <- gam[,Q1_index]
  Q3 <- gam[,Q3_index]
  Q1_psi <- sqrt(gradient(Q1,1/(M-1)))
  Q3_psi <- sqrt(gradient(Q3,1/(M-1)))

  # compute phase whiskers
  IQR <- dx[Q1_index] + dx[Q3_index]
  v1 <- v[,Q1_index]
  v3 <- v[,Q3_index]
  upper_v <- v3 + kp * IQR * v3 / sqrt(trapz(time,v3*v3))
  lower_v <- v1 + kp * IQR * v1 / sqrt(trapz(time,v1*v1))
  upper_psi <- exp_map(psi_median, upper_v)
  lower_psi <- exp_map(psi_median, lower_v)
  upper <- cumtrapz(time,upper_psi*upper_psi)
  lower <- cumtrapz(time,lower_psi*lower_psi)

  upper_dis <- sqrt(trapz(time,(upper_v)^2))
  lower_dis <- sqrt(trapz(time,(lower_v)^2))
  whisker_dis <- max(c(upper_dis,lower_dis))

  # indentify phase outliers
  outlier_index <- c()
  for (i in 1:N){
    if (dx[dx_ordering[N+1-i]]> whisker_dis){
      outlier_index <- c(outlier_index, dx_ordering[N+1-i])
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
    distance_to_upper[j] = sqrt(trapz(time,(upper_v-v[,j])^2))
    distance_to_lower[j] = sqrt(trapz(time,(lower_v-v[,j])^2))
  }
  max_index <- which.min(distance_to_upper)
  min_index <- which.min(distance_to_lower)
  min_psi <- psi[,min_index]
  max_psi <- psi[,max_index]
  minn <- gam[,min_index]
  maxx <- gam[,max_index]

  if (showplot){
    plot(time, median_x, col="black",xlab="Time",main="Phase Boxplot", type="l", ylim=c(0, 1))
    lines(time, Q1, col="blue")
    lines(time, Q3, col="green")
    lines(time, minn, col="red")
    lines(time, maxx, col="magenta")

    s <- seq(0,1,length.out=100)
    Fs2 <- matrix(0,length(time), 397)
    Fs2[,1] <- (1-s[1]) * (minn-time) + s[1] * (Q1-time)
    for (j in 2:100){
      Fs2[,j] <- (1-s[j]) * (minn-time) + s[j] * (Q1-time)
      Fs2[,99+j] <- (1-s[j]) * (Q1-time) + s[j] * (median_x-time)
      Fs2[,198+j] <- (1-s[j]) * (median_x-time) + s[j] * (Q3-time)
      Fs2[,297+j] <- (1-s[j]) * (Q3-time) + s[j] * (maxx-time)
    }
    d1<-sqrt(trapz(time,(min_psi*Q1_psi)^2))
    dl<-sqrt(trapz(time,(Q1_psi*psi_median)^2))
    d3<-sqrt(trapz(time,(psi_median*Q3_psi)^2))
    du<-sqrt(trapz(time,(Q3_psi*max_psi)^2))
    part1<-seq(-d1-dl,-d1,length.out=100)
    part2<-seq(-d1,0,length.out=100)
    part3<-seq(0,d3,length.out=100)
    part4<-seq(d3,d3+du,length.out=100)
    allparts<-c(part1,part2[2:100],part3[2:100],part4[2:100])

    image(time, allparts, Fs2, main="Surface Plot", ylab="", col=viridis(128))

  }

  return(list(median_x=median_x,Q1=Q1,Q3=Q3,minn=minn,maxx=maxx,outlier_index=outlier_index))
}
