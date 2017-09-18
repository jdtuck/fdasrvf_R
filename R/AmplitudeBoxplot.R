#' Amplitude Boxplot
#'
#' This function constructs the amplitude boxplot
#'
#' @param warp_median fdawarp object from \link{time_warping} of aligned data using the median
#' @param alpha quantile value (default=.05, i.e., 95\%)
#' @param ka scalar for outlier cutoff (default=1)
#' @param showplot shows plots of functions (default = T)
#' @return Returns a ampbox object containing \item{median_y}{median function}
#' \item{Q1}{First quartile}
#' \item{Q3}{Second quartile}
#' \item{Q1a}{First quantile based on alpha}
#' \item{Q3a}{Second quantile based on alpha}
#' \item{minn}{minimum extreme function}
#' \item{maxx}{maximum extreme function}
#' \item{outlier_index}{indexes of outlier functions}
#' \item{fmedian}{median function}
#' @keywords srvf alignment boxplot
#' @references Xie, W., S. Kurtek, K. Bharath, and Y. Sun  (2016). "A Geometric Approach to Visualization
#' of Variability in Functional Data." Journal of the American Statistical Association in press: 1-34.
#' @export
#' @examples
#' data("simu_warp_median")
#' out <- AmplitudeBoxplot(simu_warp_median, showplot=FALSE)
AmplitudeBoxplot <- function(warp_median, alpha=.05, ka=1, showplot=TRUE){

  fn <- warp_median$fn
  median_y <- warp_median$fmean
  qn <- warp_median$qn
  qmedian <- warp_median$mqn
  time <- warp_median$time
  if (warp_median$method != 'median'){
      stop('need aligned to median, please rerun time_warping with method="median"')
  }

  if (warp_median$rsamps){
    fn <- warp_median$fs
    qn <- warp_median$qs
  }

  M <- nrow(fn)
  N <- ncol(fn)
  lambda <- 0.5

  # translation
  translation <- rep(0,N)
  for (i in 1:N){
    translation[i] <- trapz(time, fn[,i]/(time[M]-time[1]))
  }

  # compute amplitude distances
  dy <- rep(0, N)
  for (i in 1:N){
    dy[i] <- sqrt(trapz(time, (qmedian-qn[,i])^2))
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

  # identify amplitude quantile
  dy_ordering <- sort(dy, index.return = T)$ix
  CR_alpha <- dy_ordering[1:round(N*(1-alpha))]  # (1-alpha)% central region
  m <- max(dy[CR_alpha])  # maximal amplitude distance with (1-alpha)% central region
  angle <- matrix(0, length(CR_alpha), length(CR_alpha))
  energy <- matrix(0, length(CR_alpha), length(CR_alpha))
  for (i in 1:(length(CR_alpha)-1)){
    for (j in (i+1):length(CR_alpha)){
      q1 <- qn[,CR_alpha[i]] - qmedian
      q3 <- qn[,CR_alpha[j]] - qmedian
      q1 <- q1/sqrt(trapz(time,q1*q1))
      q3 <- q3/sqrt(trapz(time,q3*q3))
      angle[i,j] <- trapz(time,q1*q3)
      energy[i,j] <- (1-lambda) * (dy[CR_alpha[i]]/m + dy[CR_alpha[j]]/m) - lambda * (angle[i,j]+1)
    }
  }
  maxloc <- which(energy == max(energy), arr.ind = TRUE)

  Q1a_index <- CR_alpha[maxloc[1,1]]
  Q3a_index <- CR_alpha[maxloc[1,2]]
  Q1a_q <- qn[,Q1a_index]
  Q3a_q <- qn[,Q3a_index]
  Q1a <- fn[,Q1a_index]
  Q3a <- fn[,Q3a_index]

  # compute amplitude whiskers
  IQR <- dy[Q1_index] + dy[Q3_index]
  v1 <- Q1_q - qmedian
  v3 <- Q3_q - qmedian
  upper_q <- Q3_q + ka * IQR * v3 / sqrt(trapz(time,v3*v3))
  lower_q <- Q1_q + ka * IQR * v1 / sqrt(trapz(time,v1*v1))
  upper <- cumtrapz(time,upper_q*abs(upper_q))
  lower <- cumtrapz(time,lower_q*abs(lower_q))

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
  max_index <- which.min(distance_to_upper)
  min_index <- which.min(distance_to_lower)
  min_q <- qn[,min_index]
  max_q <- qn[,max_index]
  minn <- fn[,min_index]
  maxx <- fn[,max_index]


  out <- list(median_y=median_y,Q1=Q1,Q3=Q3,Q1a=Q1a,Q3a=Q3a,minn=minn,maxx=maxx,
              outlier_index=outlier_index)
  out$fmedian <- median_y
  out$time <- time
  out$qmedian <- qmedian
  out$Q1_q <- Q1_q
  out$Q1a_q <- Q1a_q
  out$Q3_q <- Q3_q
  out$Q3a_q <- Q3a_q
  out$min_q <- min_q
  out$max_q <- max_q
  class(out) <- 'ampbox'

  if (showplot){
    plot(out)
  }

  return(out)
}
