\name{time_warping}
\alias{time_warping}
\title{Group-wise function alignment}
\usage{
  time_warping(f, time, lambda = 0, showplot = TRUE,
    smooth_data = FALSE, sparam = 25, parallel = FALSE,
    cores = 8)
}
\arguments{
  \item{f}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} functions
  with \eqn{N} samples}

  \item{time}{vector of size \eqn{N} describing the sample
  points}

  \item{lambda}{controls the elasticity (default = 0)}

  \item{showplot}{shows plots of functions (default = T)}

  \item{smooth_data}{smooth data using box filter (default
  = F)}

  \item{sparam}{number of times to apply box filter
  (default = 25)}

  \item{parallel}{enable parallel mode using
  \code{\link{foreach}} and \code{\link{doMC}} pacakge}

  \item{cores}{set number of cores to use with
  \code{\link{doMC}} (default = 8)}
}
\value{
  Returns a list containing \item{f0}{original functions}
  \item{fn}{aligned functions} \item{qn}{aligned srvfs}
  \item{q0}{original srvfs} \item{fmean}{function mean}
  \item{mqn}{srvf mean} \item{gam}{warping functions}
}
\description{
  This function aligns a collection of functions using the
  elastic square-root velocity (srvf) framework.
}
\examples{
data("simu_data")
out = time_warping(f,time)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment}
\keyword{srvf}

