#' Gaussian model of functional data using joint Model
#'
#' This function models the functional data using a Gaussian model extracted from
#' the principal components of the srvfs using the joint model
#'
#' @param warp_data fdawarp object from [time_warping] of aligned data
#' @param n number of random samples (n = 1)
#' @param no number of principal components (n=4)
#' @return Returns a fdawarp object containing \item{fs}{random aligned samples}
#' \item{gams}{random warping function samples}
#' \item{ft}{random function samples}
#' \item{qs}{random srvf samples}
#' @keywords pca
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references Jung, S. L. a. S. (2016). "Combined Analysis of Amplitude and Phase Variations in Functional Data."
#'        	arXiv:1603.01775 [stat.ME].
#' @export
#' @examples
#' out1 <- joint_gauss_model(simu_warp, n = 10)
joint_gauss_model <- function(warp_data, n=1, no=5){
  fn <- warp_data$fn
  time <- warp_data$time
  qn <- warp_data$qn
  gam <- warp_data$gam
  # Perform PCA -------------------------------------------------------------
  M <- nrow(fn)
  jfpca <- jointFPCA(warp_data, no, showplot = F)
  s <- jfpca$latent
  U <- jfpca$U
  C <- jfpca$C
  mu_psi <- jfpca$mu_psi

  mq_new <- rowMeans(qn)
  id <- round(length(time)/2)
  m_new <- sign(fn[id,])*sqrt(abs(fn[id,]))  # scaled version
  mqn <- c(mq_new,mean(m_new))

  # Generate Random Samples -------------------------------------------------
  if (length(s)>1){
    vals <- rmvnorm(n, sigma = diag(s))
  } else {
    vals <- rnorm(n, sd=s)
  }

  tmp <- U%*%t(vals)
  qhat <- matrix(mqn,M+1,n) + tmp[1:(M+1),]
  tmp <- U%*%t(vals/C)
  vechat <- tmp[(M+2):nrow(tmp),]
  psihat <- apply(vechat,2,function(x,mu_psi){exp_map(mu_psi,x)},mu_psi)
  gamhat <- apply(psihat,2,function(x,time){gam_mu = cumtrapz(time, x*x)
                  gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))},
                  seq(0,1,length.out=M))

  ft <- array(0,dim=c(M,n))
  fhat <- array(0,dim=c(M,n))
  for (ii in 1:n){
    fhat[,ii] <- cumtrapzmid(time,qhat[1:M,ii]*abs(qhat[1:M,ii]),sign(qhat[M+1,ii])*(qhat[M+1,ii]^2), id)
    ft[,ii] <- warp_f_gamma(fhat[,ii], seq(0,1,length.out=M), gamhat[,ii])
  }

  warp_data$fs <- fhat
  warp_data$gams <- gamhat
  warp_data$ft <- ft
  warp_data$qs <- qhat[1:M,]
  warp_data$rsamps <- T

  return(warp_data)
}
