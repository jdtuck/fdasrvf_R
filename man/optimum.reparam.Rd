\name{optimum.reparam}
\alias{optimum.reparam}
\title{Align two functions}
\usage{
  optimum.reparam(Q1, T1, Q2, T2)
}
\arguments{
  \item{Q1}{srvf of function 1}

  \item{T1}{sample points of function 1}

  \item{Q2}{srvf of function 2}

  \item{T2}{sample points of function 2}
}
\value{
  gam warping function
}
\description{
  This function aligns two srvf functions using Dynamic
  Programming
}
\examples{
data("simu_data")
q = f_to_srvf(simu_data$f,simu_data$time)
gam = optimum.reparam(q[,1],simu_data$time,q[,2],simu_data$time)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment,}
\keyword{pca}
\keyword{srvf}

