\name{align_fPCA}
\alias{align_fPCA}
\title{Group-wise function alignment and PCA Extractions}
\usage{
  align_fPCA(f, time, num_comp = 3, showplot = T,
    smooth_data = FALSE, sparam = 25, parallel = FALSE,
    cores = 8)
}
\arguments{
  \item{f}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} functions
  with \eqn{N} samples}

  \item{time}{vector of size \eqn{N} describing the sample
  points}

  \item{num_comp}{number of principal components to extract
  (default = 3)}

  \item{showplot}{shows plots of functions (default = T)}

  \item{smooth_data}{smooth data using box filter (default
  = F)}

  \item{sparam}{number of times to apply box filter
  (default = 25)}

  \item{parallel}{enable parallel mode using
  \code{\link{foreach}} and \code{doMC} pacakge if on linux
  or mac or \code{doSNOW} if on windows}

  \item{cores}{set number of cores to use with \code{doMC}
  or \code{doSNOW}(default = 2)}
}
\value{
  Returns a list containing \item{f0}{original functions}
  \item{fn}{aligned functions} \item{qn}{aligned srvfs}
  \item{q0}{original srvfs} \item{mqn}{srvf mean}
  \item{gam}{warping functions} \item{Dx}{cost function}
  \item{vfpca}{list containing} \item{q_pca}{srvf principal
  directions} \item{f_pca}{f principal directions}
  \item{latent}{latent values} \item{coef}{coefficients}
  \item{U}{eigenvectors}
}
\description{
  This function aligns a collection of functions while
  extracting pincipal components.
}
\examples{
data("simu_data")
out = align_fPCA(simu_data$f,simu_data$time)
}
\references{
  Tucker, J. D., Wu, W., Srivastava, A., Generative Models
  for Function Data using Phase and Amplitude Separation,
  Computational Statistics and Data Analysis (2012),
  10.1016/j.csda.2012.12.001.
}
\keyword{alignment,}
\keyword{pca}
\keyword{srvf}

