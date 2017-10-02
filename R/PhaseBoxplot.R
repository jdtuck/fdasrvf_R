#' Phase Boxplot
#'
#' This function constructs the amplitude boxplot
#'
#' @param warp_median fdawarp object from \link{time_warping} of aligned data using the median
#' @param alpha quantile value (default=.05, i.e., 95\%)
#' @param kp scalar for outlier cutoff (default=1)
#' @param showplot shows plots of functions (default = T)
#' @return Returns a phbox object containing \item{median_x}{median warping function}
#' \item{Q1}{First quartile}
#' \item{Q3}{Second quartile}
#' \item{Q1a}{First quantile based on alpha}
#' \item{Q3a}{Second quantile based on alpha}
#' \item{minn}{minimum extreme function}
#' \item{maxx}{maximum extreme function}
#' \item{outlier_index}{indexes of outlier functions}
#' @keywords srvf alignment boxplot
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A Geometric Approach to Visualization
#' of Variability in Functional Data." Journal of the American Statistical Association in press: 1-34.
#' @export
#' @examples
#' data("simu_warp_median")
#' out <- PhaseBoxplot(simu_warp_median, showplot=FALSE)
PhaseBoxplot <- function(warp_median, alpha=.05, kp=1, showplot=TRUE){
  gam <- warp_median$gam
  if (warp_median$method != 'median'){
      stop('need aligned to median, please rerun time_warping with method="median"')
  }

  if (warp_median$rsamps){
    gam <- warp_median$gams
  }

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
    dx[i] <- sqrt(trapz(time, v[,i]^2))
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

  # identify phase quantiles
  dx_ordering <- sort(dx, index.return = T)$ix
  CR_alpha <- dx_ordering[1:round(N*(1-alpha))]  # (1-alpha)% central region
  m <- max(dx[CR_alpha])  # maximal phase distance with (1-alpha)% central region
  angle <- matrix(0, length(CR_alpha), length(CR_alpha))
  energy <- matrix(0, length(CR_alpha), length(CR_alpha))
  for (i in 1:(length(CR_alpha)-1)){
    for (j in (i+1):length(CR_alpha)){
      q1 <- v[,CR_alpha[i]]
      q3 <- v[,CR_alpha[j]]
      q1 <- q1/sqrt(trapz(time,q1*q1))
      q3 <- q3/sqrt(trapz(time,q3*q3))
      angle[i,j] <- trapz(time,q1*q3)
      energy[i,j] <- (1-lambda) * (dx[CR_alpha[i]]/m + dx[CR_alpha[j]]/m) - lambda * (angle[i,j]+1)
    }
  }
  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1a_index <- CR_alpha[maxloc[1,1]]
  Q3a_index <- CR_alpha[maxloc[1,2]]
  Q1a <- gam[,Q1a_index]
  Q3a <- gam[,Q3a_index]
  Q1a_psi <- sqrt(gradient(Q1a,1/(M-1)))
  Q3a_psi <- sqrt(gradient(Q3a,1/(M-1)))

  # check quartile and quantile going same direction
  tst <- trapz(time,v[,Q1a_index]*v[,Q1_index])
  if (tst < 0){
    Q1a <- gam[,Q3a_index]
    Q3a <- gam[,Q1a_index]
  }

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

  out <- list(median_x=median_x,Q1=Q1,Q3=Q3,Q1a=Q1a,Q3a=Q3a,minn=minn,maxx=maxx,
              outlier_index=outlier_index)
  out$time <- time
  out$psi_median <- psi_median
  out$Q1_psi <- Q1_psi
  out$Q1a_psi <- Q1a_psi
  out$Q3_psi <- Q3_psi
  out$Q3a_psi <- Q3a_psi
  out$max_psi <- max_psi
  out$min_psi <- min_psi
  class(out) <- 'phbox'

  if (showplot){
    plot(out)

  }

  return(out)
}
