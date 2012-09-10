\name{f_to_srvf}
\alias{f_to_srvf}
\title{Convert to SRVF}
\usage{
  f_to_srvf(f, time)
}
\arguments{
  \item{f}{matrix of functions}

  \item{time}{time}
}
\value{
  q matrix of srvfs
}
\description{
  This function converts functions to srvf
}
\examples{
data("simu_data")
q = f_to_srvf(f,time)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment}
\keyword{srvf}

