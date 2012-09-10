\name{SqrtMean}
\alias{SqrtMean}
\title{SRVF transform of warping functions}
\usage{
  SqrtMean(gam)
}
\arguments{
  \item{gam}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} warping
  functions}
}
\value{
  Returns a list containing \item{mu}{mean function}
  \item{psi}{srvf of warping functions} \item{vec}{shooting
  vectors}
}
\description{
  This function calculates the srvf of warping functions
  with corresponding shooting vectors
}
\examples{
data("simu_data")
gam = time_warping(simu_data$f,simu_data$time)$gam
out = SqrtMean(gam)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment}
\keyword{srvf}

