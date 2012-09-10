\name{gauss_model}
\alias{gauss_model}
\title{Gaussian model of functional data}
\usage{
  gauss_model(fn, time, qn, gam, n = 1,
    sort_samples = FALSE)
}
\arguments{
  \item{fn}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} aligned
  functions with \eqn{N} samples}

  \item{time}{vector of size \eqn{N} describing the sample
  points}

  \item{qn}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} aligned
  srvfs}

  \item{gam}{warping functions}

  \item{n}{number of random samples (n = 1)}

  \item{sort_samples}{sort samples (default = F)}
}
\value{
  Returns a list containing \item{fs}{random aligned
  samples} \item{gams}{random warping function samples}
}
\description{
  This function models the functional data using a Gaussian
  model extracted from the principal components of the
  srvfs
}
\examples{
data("simu_data")
out = time_warping(f,time)
out1 = gauss_model(out$fn,time,out$qn,out$gam,n = 10)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{pca}

