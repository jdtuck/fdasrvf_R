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
data("simu_data")
fo = smooth.data(simu_data$f,25)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].

  Tucker, J. D., Wu, W., Srivastava, A., Dec 2012.
  Generative Models for Function Data using Phase and
  Amplitude Separation, accepted to Computational
  Statistics and Data Analysis.
}
\keyword{alignment,}
\keyword{pca}
\keyword{srvf}

