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
q = f_to_srvf(simu_data$f,simu_data$time)
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
\keyword{alignment}
\keyword{srvf}

