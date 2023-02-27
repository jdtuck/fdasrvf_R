#' Tolerance Bound Calculation using Bootstrap Sampling
#'
#' This function computes tolerance bounds for functional data containing
#' phase and amplitude variation using bootstrap sampling
#'
#' @param f matrix of functions
#' @param time vector describing time sampling
#' @param a confidence level of tolerance bound (default = 0.05)
#' @param p coverage level of tolerance bound (default = 0.99)
#' @param B number of bootstrap samples (default = 500)
#' @param no number of principal components (default = 5)
#' @param Nsamp number of functions per bootstrap (default = 100)
#' @param parallel enable parallel processing (default = T)
#' @return Returns a list containing \item{amp}{amplitude tolerance bounds}
#' \item{ph}{phase tolerance bounds}
#' @keywords tolerance bootstrap
#' @concept bounds
#' @references J. D. Tucker, J. R. Lewis, C. King, and S. Kurtek, “A Geometric
#'   Approach for Computing Tolerance Bounds for Elastic Functional Data,”
#'   Journal of Applied Statistics, 10.1080/02664763.2019.1645818, 2019.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative Models for
#'   Function Data using Phase and Amplitude Separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references Jung, S. L. a. S. (2016). "Combined Analysis of Amplitude and
#'   Phase Variations in Functional Data." arXiv:1603.01775.
#' @export
#' @examples
#' \dontrun{
#'   out1 <- bootTB(simu_data$f, simu_data$time)
#' }
bootTB <- function(f, time, a=.05, p=.99, B=500, no = 5, Nsamp=100, parallel=T){

  M <- nrow(f)
  N <- ncol(f)
  # Align Data --------------------------------------------------------------
  out.med <- time_warping(
    f, time,
    centroid_type = "median",
    parallel = parallel,
    show_plot = FALSE
  )

  if (parallel){
    cores <- max(parallel::detectCores() - 1, 1)
    if (cores > 40){
      cores = 32
    }
    cl <- parallel::makeCluster(cores, outfile = "")
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  # Caclculate CI -----------------------------------------------------------
  # a% tolerance bound with p% coverage
  cat("Bootstrap Sampling\n")
  k = 1
  pb <- utils::txtProgressBar(0, B, style = 3)
  outfor <- foreach::foreach(k=1:B, .combine=cbind, .packages=c('fdasrvf','mvtnorm')) %dopar% {
    samples <- joint_gauss_model(out.med, Nsamp, no)
    amp <- AmplitudeBoxplot(samples, alpha=1-p, showplot=F)
    ph <- PhaseBoxplot(samples, alpha=1-p, showplot=F)
    utils::setTxtProgressBar(pb, k)
    list(amp$Q3a,amp$Q1a,ph$Q3a,ph$Q1a)
  }
  bootLwr.amp <- unlist(outfor[1,])
  dim(bootLwr.amp)=c(M,B)
  bootUpr.amp <- unlist(outfor[2,])
  dim(bootUpr.amp)=c(M,B)
  bootLwr.ph <- unlist(outfor[3,])
  dim(bootLwr.ph)=c(M,B)
  bootUpr.ph <- unlist(outfor[4,])
  dim(bootUpr.ph)=c(M,B)

  boot.amp <- cbind(bootLwr.amp, bootUpr.amp)
  boot.ph <- cbind(bootLwr.ph, bootUpr.ph)
  boot.amp.q <- f_to_srvf(boot.amp, time)
  boot.out <- out.med
  boot.out$fn <- boot.amp
  boot.out$qn <- boot.amp.q
  boot.out$gam <- boot.ph

  if (parallel) parallel::stopCluster(cl)

  # Tolerance Bounds --------------------------------------------------------
  amp <- AmplitudeBoxplot(boot.out, alpha=a, showplot=F)
  ph <- PhaseBoxplot(boot.out, alpha=a, showplot=F)

  return(list(amp=amp,ph=ph,align=out.med))

}

#' Tolerance Bound Calculation using Elastic Functional PCA
#'
#' This function computes tolerance bounds for functional data containing
#' phase and amplitude variation using principal component analysis
#'
#' @param f matrix of functions
#' @param time vector describing time sampling
#' @param m number of principal components (default = 4)
#' @param B number of monte carlo iterations
#' @param a confidence level of tolerance bound (default = 0.05)
#' @param p coverage level of tolerance bound (default = 0.99)
#' @return Returns a list containing \item{pca}{pca output}
#' \item{tol}{tolerance factor}
#' @concept pca tolerance bounds
#' @references J. D. Tucker, J. R. Lewis, C. King, and S. Kurtek, “A Geometric
#'   Approach for Computing Tolerance Bounds for Elastic Functional Data,”
#'   Journal of Applied Statistics, 10.1080/02664763.2019.1645818, 2019.
#' @references Tucker, J. D., Wu, W., Srivastava, A., Generative Models for
#'   Function Data using Phase and Amplitude Separation, Computational
#'   Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @references Jung, S. L. a. S. (2016). "Combined Analysis of Amplitude and
#'   Phase Variations in Functional Data." arXiv:1603.01775.
#' @export
#' @examples
#' \dontrun{
#'   out1 <- pcaTB(simu_data$f, simu_data$time)
#' }
pcaTB <- function(f, time, m = 4, B = 100000, a = 0.05, p = 0.99){
  # Align Data --------------------------------------------------------------
  out <- time_warping(f, time, parallel = TRUE, show_plot = FALSE)

  # Calculate PCA -----------------------------------------------------------
  out.pca <- jointFPCA(out, m, showplot=F)

  # Calcualte TB
  # Krishnamoorthy, K. and Mondal, S. (2006), Improved Tolerance Factors for Multivariate Normal
  # Distributions, Communications in Statistics - Simulation and Computation, 35, 461–478.
  tol <- tolerance::mvtol.region(x = out.pca$coef, alpha = a, P = p, B = B)

  return(list(pca=out.pca, align=out, tol=tol))
}
