#' K-Means Clustering and Alignment
#'
#' This function clusters functions and aligns using the elastic square-root
#' slope (srsf) framework.
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param K number of clusters
#' @param seeds indexes of cluster center functions (default = NULL)
#' @param lambda controls the elasticity (default = 0)
#' @param showplot shows plots of functions (default = T)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using \code{\link{foreach}} and
#'   \code{doParallel} pacakge (default=F)
#' @param alignment wether to perform alignment (default = T)
#' @param omethod optimization method (DP,DP2,RBFGS)
#' @param MaxItr maximum number of iterations
#' @param thresh cost function threshold
#' @return Returns a fdakma object containing \item{f0}{original functions}
#' \item{fn}{aligned functions - matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples which is a list for each cluster}
#' \item{qn}{aligned SRSFs - similar structure to fn}
#' \item{q0}{original SRSFs}
#' \item{labels}{cluster labels}
#' \item{templates}{cluster center functions}
#' \item{templates.q}{cluster center SRSFs}
#' \item{gam}{warping functions - similar structure to fn}
#' \item{qun}{Cost Function Value}
#' @keywords srsf alignment clustering
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references Sangalli, L. M., et al. (2010). "k-mean alignment for curve clustering."
#'  Computational Statistics & Data Analysis 54(5): 1219-1233.
#' @export
#' @examples
#' data("growth_vel")
#' # use more iterations for accuracy
#' out <- kmeans_align(growth_vel$f,growth_vel$time, K=2, MaxItr=1)
kmeans_align <- function(f, time, K, seeds=NULL, lambda = 0, showplot = TRUE,
                   smooth_data = FALSE, sparam = 25, parallel = FALSE,
                   alignment = TRUE, omethod = "DP", MaxItr = 50, thresh = 0.01){
  # Initialize --------------------------------------------------------------
  w <- 0.0
  k <- 1
  if (parallel){
    cores <- detectCores()
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  } else
  {
    registerDoSEQ()
  }

  M <- nrow(f)
  N <- ncol(f)
  if (is.null(seeds)){
    template.ind <- sample(1:N,K)
  } else {
    template.ind <- seeds
  }
  templates <- matrix(0,M,K)
  for (i in 1:K){
    templates[,i] <- f[,template.ind[i]]
  }
  cluster.id <- rep(0,N)
  qun <- rep(0, MaxItr)

  # Convert to SRSF
  if (smooth_data){
    f <- smooth.data(f,sparam)
  }

  if (showplot){
    matplot(time,f,type="l")
    title(main="Original data")
  }

  q <- f_to_srvf(f, time)
  templates.q <- matrix(0,M,K)
  for (i in 1:K){
    templates.q[,i] <- q[,template.ind[i]]
  }

  for (itr in 1:MaxItr){
    cat(sprintf("updating step: r=%d\n", itr))
    # Assignment and Alignment ------------------------------------------------
    gam <- list()
    Dy <- matrix(0,K,N)
    qn <- list()
    fn <- list()
    for (i in 1:K){
      outfor<-foreach(k = 1:N, .combine=cbind, .packages="fdasrvf") %dopar% {
        if (alignment){
          gam_tmp <- optimum.reparam(templates.q[,i],time,q[,k],time,lambda,omethod,w,templates[1,i],f[1,k])
        } else {
          gam_tmp <- seq(0,1,length.out=M)
        }
        fw <- approx(time,f[,k],xout=(time[length(time)]-time[1])*gam_tmp + time[1])$y
        qw <- f_to_srvf(fw,time)
        dist <- sqrt(sum(trapz(time, (qw-templates.q[,i])^2)))
        list(gam_tmp,fw,qw,dist)
      }
      gamt <- unlist(outfor[1,]);
      dim(gamt) <- c(M,N)
      gam[[i]] <- gamt
      f_temp <- unlist(outfor[2,]);
      dim(f_temp) <- c(M,N)
      q_temp <- unlist(outfor[3,]);
      dim(q_temp) <- c(M,N)
      qn[[i]] <- q_temp
      fn[[i]] <- f_temp
      dtil <- unlist(outfor[4,]);
      dim(dtil) <- c(1,N)
      Dy[i,] <- dtil
    }
    cluster.id <- apply(Dy,2,which.min)


    # Normalization -----------------------------------------------------------
    for (i in 1:K){
      id <- which(cluster.id==i)
      ftmp <- fn[[i]][,id]
      gamtmp <- gam[[i]][,id]
      gamI <- SqrtMeanInverse(gamtmp)
      N1 <- length(id)
      outfor<-foreach(k = 1:N1, .combine=cbind, .packages="fdasrvf") %dopar% {
        fw <- approx(time,ftmp[,k],xout=(time[length(time)]-time[1])*gamI +
                       time[1])$y
        qw <- f_to_srvf(fw,time)
        gamt1 <- approx(time,gamtmp[,k],xout=(time[length(time)]-time[1])*gamI +
                        time[1])$y
        list(gamt1,fw,qw)
      }
      gamt <- unlist(outfor[1,]);
      dim(gamt) <- c(M,N1)
      gam[[i]][,id] <- gamt
      f_temp <- unlist(outfor[2,]);
      dim(f_temp) <- c(M,N1)
      q_temp <- unlist(outfor[3,]);
      dim(q_temp) <- c(M,N1)
      qn[[i]][,id] <- q_temp
      fn[[i]][,id] <- f_temp
    }


    # Template Identification -------------------------------------------------
    qun.t <- rep(0,K)
    for (i in 1:K){
      id <- which(cluster.id==i)
      old.templates.q <- templates.q
      templates.q[,i] <- rowMeans(qn[[i]][,id])
      templates[,i] <- rowMeans(fn[[i]][,id])

      qun.t[i] <- pvecnorm(templates.q[,i]-old.templates.q[,i],2)/pvecnorm(old.templates.q[,i],2)
    }
    qun[itr] <- mean(qun.t)

    if (qun[itr] < thresh)
      break
  }


  # Output ------------------------------------------------------------------
  ftmp <- qtmp <- gamtmp <- list()
  for (i in 1:K){
    id <- which(cluster.id==i)
    ftmp[[i]] <- fn[[i]][,id]
    qtmp[[i]] <- qn[[i]][,id]
    gamtmp[[i]] <- gam[[i]][,id]
  }

  out <- list(f0=f,q0=q,time=time,fn=ftmp,qn=qtmp,gam=gamtmp,labels=cluster.id,
              templates=templates,templates.q=templates.q,lambda=lambda,omethod=omethod,
              qun=qun[1:itr])

  class(out) <- 'fdakma'

  if (showplot){
    plot(out)
  }

  if (parallel){
    stopCluster(cl)
  }

  return(out)
}
