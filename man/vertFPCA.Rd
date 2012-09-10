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
out = time_warping(f,time)
vfpca = vertFPCA(out$fn,time,out$qn,no = 3)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment}
\keyword{srvf}

