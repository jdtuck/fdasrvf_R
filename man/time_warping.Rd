\name{time_warping}
\alias{time_warping}
\title{Group-wise function alignment}
\usage{
  time_warping(f, time, lambda = 0, method = "mean",
    showplot = TRUE, smooth_data = FALSE, sparam = 25,
    parallel = FALSE, cores = 2)
}
\arguments{
  \item{f}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} functions
  with \eqn{N} samples}

  \item{time}{vector of size \eqn{N} describing the sample
  points}

  \item{lambda}{controls the elasticity (default = 0)}

  \item{method}{warp and caluclate to Karcher Mean or
  Median (options = "mean" or "median", default = "mean")}

  \item{showplot}{shows plots of functions (default = T)}

  \item{smooth_data}{smooth data using box filter (default
  = F)}

  \item{sparam}{number of times to apply box filter
  (default = 25)}

  \item{parallel}{enable parallel mode using
  \code{\link{foreach}} and \code{doMC} pacakge if on linux
  or mac or \code{doSNOW} if on windows}

  \item{cores}{set number of cores to use with \code{doMC}
  or \code{doSNOW}(default = 2)}
}
\value{
  Returns a list containing \item{f0}{original functions}
  \item{fn}{aligned functions} \item{qn}{aligned srvfs}
  \item{q0}{original srvfs} \item{fmean}{function mean or
  median} \item{mqn}{srvf mean or median}
  \item{gam}{warping functions} \item{orig.var}{Original
  Variance of Functions} \item{amp.var}{Amplitude Variance}
  \item{phase.var}{Phase Variance}
}
\description{
  This function aligns a collection of functions using the
  elastic square-root velocity (srvf) framework.
}
\examples{
data("simu_data")
out = time_warping(simu_data$f,simu_data$time)
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

