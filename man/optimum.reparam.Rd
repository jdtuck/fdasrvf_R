\name{optimum.reparam}
\alias{optimum.reparam}
\title{Align two functions}
\usage{
  optimum.reparam(Q1, T1, Q2, T2, lambda = 0)
}
\arguments{
  \item{Q1}{srvf of function 1}

  \item{T1}{sample points of function 1}

  \item{Q2}{srvf of function 2}

  \item{T2}{sample points of function 2}

  \item{lambda}{controls amount of warping (default = 0)}
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

  Tucker, J. D., Wu, W., Srivastava, A., Dec 2012.
  Generative Models for Function Data using Phase and
  Amplitude Separation, accepted to Computational
  Statistics and Data Analysis.
}
\keyword{alignment,}
\keyword{pca}
\keyword{srvf}

