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
  \code{\link{foreach}} and \code{\link{doMC}} pacakge}

  \item{cores}{set number of cores to use with
  \code{\link{doMC}} (default = 8)}
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
out = align_fPCA(f,time)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment,}
\keyword{pca}
\keyword{srvf}

