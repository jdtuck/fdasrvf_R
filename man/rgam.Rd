\name{rgam}
\alias{rgam}
\title{Random Warping}
\usage{
  rgam(N, sigma, num)
}
\arguments{
  \item{N}{length of warping function}

  \item{sigma}{variance of warping functions}

  \item{num}{number of warping functions}
}
\value{
  gam warping functions
}
\description{
  Generates random warping functions
}
\examples{
gam = rgam(N=101, sigma=.01, num=35)
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
\keyword{diffeomorphism,}
\keyword{function}
\keyword{warping}

