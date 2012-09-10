\name{smooth.data}
\alias{smooth.data}
\title{Smooth Functions}
\usage{
  smooth.data(f, sparam)
}
\arguments{
  \item{f}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} functions
  with \eqn{N} samples}

  \item{sparam}{number of times to run box filter}
}
\value{
  fo smoothed functions
}
\description{
  This function smooths functions using standard box filter
}
\examples{
fo = smooth.data(f,25)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment,}
\keyword{pca}
\keyword{srvf}

