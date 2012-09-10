\name{vertFPCA}
\alias{vertFPCA}
\title{Vertical Functional Principal Component Analysis}
\usage{
  vertFPCA(fn, time, qn, no, showplot = TRUE)
}
\arguments{
  \item{fn}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} aligned
  functions with \eqn{N} samples}

  \item{time}{vector of size \eqn{N} describing the sample
  points}

  \item{qn}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} of
  aligned srvfs}

  \item{no}{number of prinicpal components to extract}

  \item{showplot}{show plots of prinipal directions
  (default = T)}
}
\value{
  Returns a list containing \item{q_pca}{srvf principal
  directions} \item{f_pca}{f principal directions}
  \item{latent}{latent values} \item{coef}{coefficients}
  \item{U}{eigenvectors}
}
\description{
  This function calculates vertical functional principal
  component analysis on aligned data
}
\examples{
data("simu_data")
out = time_warping(simu_data$f,simu_data$time)
vfpca = vertFPCA(out$fn,simu_data$time,out$qn,no = 3)
}
\references{
  Tucker, J. D., Wu, W., Srivastava, A., May 2012.
  Generative Models for Function Data using Phase and
  Amplitude Separation, submitted to Computational
  Statistics and Data Analysis.
}
\keyword{alignment}
\keyword{srvf}

