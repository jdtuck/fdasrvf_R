% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_karcher_mean.R
\name{multivariate_karcher_mean}
\alias{multivariate_karcher_mean}
\title{Karcher Mean of Multivariate Functional Data}
\usage{
multivariate_karcher_mean(
  beta,
  mode = "O",
  alignment = TRUE,
  rotation = FALSE,
  scale = FALSE,
  lambda = 0,
  maxit = 20L,
  ms = c("mean", "median"),
  exact_medoid = FALSE,
  ncores = 1L,
  verbose = FALSE
)
}
\arguments{
\item{beta}{A numeric array of shape \eqn{L \times M \times N} specifying an
\eqn{N}-sample of \eqn{L}-dimensional functional data evaluated on a same
grid of size \eqn{M}.}

\item{mode}{A character string specifying whether the input curves should be
considered open (\code{mode == "O"}) or closed (\code{mode == "C"}). Defaults to
\code{"O"}.}

\item{alignment}{A boolean value specifying whether the curves should be
aligned before computing the distance matrix. Defaults to \code{TRUE}.}

\item{rotation}{A boolean specifying whether the distance should be made
invariant by rotation. Defaults to \code{TRUE}.}

\item{scale}{A boolean specifying whether the distance should be made
invariant by scaling. This is effectively achieved by making SRVFs having
unit length and switching to an appropriate metric on the sphere between
normalized SRVFs. Defaults to \code{TRUE}.}

\item{lambda}{A numeric value specifying the weight of a penalty term that
constraints the warping function to be not too far from the identity.
Defaults to \code{0.0}.}

\item{maxit}{An integer value specifying the maximum number of iterations.
Defaults to \code{20L}.}

\item{ms}{A character string specifying whether the Karcher mean ("mean") or
Karcher median ("median") is returned. Defaults to \code{"mean"}.}

\item{exact_medoid}{A boolean specifying whether to compute the exact medoid
from the distance matrix or as the input curve closest to the pointwise
mean. Defaults to \code{FALSE} for saving computational time.}

\item{ncores}{An integer value specifying the number of cores to use for
parallel computation. Defaults to \code{1L}. The maximum number of available
cores is determined by the \strong{parallel} package. One core is always left
out to avoid overloading the system.}

\item{verbose}{A boolean specifying whether to print the progress of the
algorithm. Defaults to \code{FALSE}.}
}
\value{
A fdacurve object with the following components:
\itemize{
\item \code{beta}: A numeric array of shape \eqn{L \times M \times N} storing the
original input data.
\item \code{q}: A numeric array of shape \eqn{L \times M \times N} storing the SRVFs
of the input data.
\item \code{betan}: A numeric array of shape \eqn{L \times M \times N} storing the
aligned, possibly optimally rotated and optimally scaled, input data.
\item \code{qn}: A numeric array of shape \eqn{L \times M \times N} storing the SRVFs
of the aligned, possibly optimally rotated and optimally scaled, input data.
\item \code{gamma}: A numeric array of shape \eqn{L \times M \times N} storing the warping
functions of the aligned, possibly optimally rotated and optimally scaled, input data.
\item \code{R}: A numeric array of shape \eqn{L \times L \times N} storing the rotation
matrices of the aligned, possibly optimally rotated and optimally scaled, input data.
\item \code{betamean}: A numeric array of shape \eqn{L \times M} storing the Karcher
mean or median of the input data.
\item \code{qmean}: A numeric array of shape \eqn{L \times M} storing the Karcher mean
or median of the SRVFs of the input data.
\item \code{type}: A character string indicating whether the Karcher mean or median
has been returned.
\item \code{E}: A numeric vector storing the energy of the Karcher mean or median at
each iteration.
\item \code{qun}: A numeric vector storing the cost function of the Karcher mean or
median at each iteration.
}
}
\description{
Calculates the Karcher mean or median of a collection of multivariate
functional data using the elastic square-root velocity (SRVF) framework.
While most of the time, the setting does not require a metric that is
invariant to rotation and scale, this can be achieved through the optional
arguments \code{rotation} and \code{scale}.
}
\examples{
out <- multivariate_karcher_mean(beta[, , 1, 1:2], maxit = 2)
# note: use more functions, small for speed
}
\references{
Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape
analysis of elastic curves in euclidean spaces. IEEE Transactions on
Pattern Analysis and Machine Intelligence, \strong{33} (7), 1415-1428.
}
\keyword{alignment}
\keyword{srvf}
