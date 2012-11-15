\name{outlier.detection}
\alias{outlier.detection}
\title{Outlier Detection}
\usage{
  outlier.detection(q, time, mq, k = 1.5)
}
\arguments{
  \item{q}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} SRVF
  functions with \eqn{N} samples}

  \item{time}{vector of size \eqn{N} describing the sample
  points}

  \item{mq}{median calcuated using
  \code{\link{time_warping}}}

  \item{k}{cutoff threshold (default = 1.5)}
}
\value{
  \item{q_outlier}{outlier functions}
}
\description{
  This function calculates outlier's using geodesic
  distances of the SRVFs from the median
}
\examples{
data("toy_data")
data("toy_warp")
q_outlier = outlier.detection(toy_warp$q0,toy_data$time,toy_warp$mqn,k=.1)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{detection}
\keyword{outlier}
\keyword{srvf}

