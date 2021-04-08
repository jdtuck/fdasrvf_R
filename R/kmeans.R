#' K-Means Clustering and Alignment
#'
#' This function clusters functions and aligns using the elastic square-root
#' slope (srsf) framework.
#'
#' @param f matrix (\eqn{N} x \eqn{M}) of \eqn{M} functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param K number of clusters
#' @param seeds indexes of cluster center functions (default = NULL)
#' @param nonempty minimum number of functions per cluster in assignment step of
#'   k-means. Set it as a positive integer to avoid the problem of empty
#'   clusters (default = 0)
#' @param lambda controls the elasticity (default = 0)
#' @param showplot shows plots of functions (default = T)
#' @param smooth_data smooth data using box filter (default = F)
#' @param sparam number of times to apply box filter (default = 25)
#' @param parallel enable parallel mode using \code{\link{foreach}} and
#'   \code{doParallel} package (default=F)
#' @param alignment whether to perform alignment (default = T)
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
#' \dontrun{
#' data("growth_vel")
#' out <- kmeans_align(growth_vel$f,growth_vel$time, K=2)
#' }
kmeans_align <- function(f, time, K, seeds=NULL, nonempty = 0, lambda = 0, 
                         showplot = TRUE, smooth_data = FALSE, sparam = 25, 
                         parallel = FALSE, alignment = TRUE, omethod = "DP", 
                         MaxItr = 50, thresh = 0.01){
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
  
  if(nonempty){
    lp.ind <- c(rbind(rep(0:(N-1), each = K), rep(N:(N+K-1), times=N)))
    lp.beg <- seq(0, N*K*2-1, 2)
    lp.rhs <- c(rep(1, N), rep(nonempty, K))
    lp.val <- rep(1,N*K*2)
    constr.ncol <- N * K
    constr.nrow <- N + K
  
    constraints.coefs <- 1*as.matrix(sparseMatrix(i = lp.ind,
                                                  p =  c(lp.beg, N*K*2), 
                                                  dims=c(constr.nrow, constr.ncol),
                                                  index1=FALSE))
    constraints.dense <- matrix(c(lp.ind+1, rep(1:(N*K), each=2), lp.val), ncol=3)
    constraints.directions <- c(rep("=", N), rep(">=", K))
    constraints.rhs <- lp.rhs
  }

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
    # Alignment ---------------------------------------------------------------
    gam <- list()
    Dy <- matrix(0,K,N)
    qn <- list()
    fn <- list()
    for (k in 1:K){
      outfor<-foreach(i = 1:N, .combine=cbind, .packages="fdasrvf") %dopar% {
        if (alignment){
          gam_tmp <- optimum.reparam(templates.q[,k],time,q[,i],time,lambda,omethod,w,templates[1,k],f[1,i])
        } else {
          gam_tmp <- seq(0,1,length.out=M)
        }
        fw <- approx(time,f[,i],xout=(time[length(time)]-time[1])*gam_tmp + time[1])$y
        qw <- f_to_srvf(fw,time)
        dist <- sqrt(sum(trapz(time, (qw-templates.q[,k])^2)))
        list(gam_tmp,fw,qw,dist)
      }
      gamt <- unlist(outfor[1,]);
      dim(gamt) <- c(M,N)
      gam[[k]] <- gamt
      f_temp <- unlist(outfor[2,]);
      dim(f_temp) <- c(M,N)
      q_temp <- unlist(outfor[3,]);
      dim(q_temp) <- c(M,N)
      qn[[k]] <- q_temp
      fn[[k]] <- f_temp
      dtil <- unlist(outfor[4,]);
      dim(dtil) <- c(1,N)
      Dy[k,] <- dtil
    }
    
    # Assignment --------------------------------------------------
    if (!nonempty) {
      cluster.id <- apply(Dy,2,which.min)
    }else{
      lpSolution <- lp(direction = "min", objective.in=c(Dy^2), 
                       const.dir = constraints.directions, const.rhs = constraints.rhs, 
                       all.bin=TRUE, dense.const=constraints.dense)
      clusterIndicator <- matrix(lpSolution$solution, nrow=K)
      if(lpSolution$status){
        stop("lpSolve failed to find feasible solution")
      }
      if(!all(clusterIndicator %in% c(0,1))){
        cat("Matrix of cluster indicators is")
        print(clusterIndicator)
        stop("Not all entries in cluster indicator matrix are binary")
      }
      cluster.id <- apply(clusterIndicator,2, which.max)
    }


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
