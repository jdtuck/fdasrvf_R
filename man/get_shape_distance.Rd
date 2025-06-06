% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refactoring.R
\name{get_shape_distance}
\alias{get_shape_distance}
\title{Computes the distance between two shapes}
\usage{
get_shape_distance(
  q1fun,
  q2fun,
  alignment = FALSE,
  rotation = FALSE,
  scale = FALSE,
  lambda = 0,
  M = 200L
)
}
\arguments{
\item{q1fun}{A function that takes a numeric vector \eqn{s} of values in
\eqn{[0, 1]} as input and returns the values of the first SRVF at \eqn{s}.}

\item{q2fun}{A function that takes a numeric vector \eqn{s} of values in
\eqn{[0, 1]} as input and returns the values of the second SRVF at \eqn{s}.}

\item{alignment}{A boolean value specifying whether the two SRVFs should be
optimally aligned before computing the distance. Defaults to \code{FALSE}.}

\item{rotation}{A boolean value specifying whether the two SRVFs should be
optimally rotated before computing the distance. Defaults to \code{FALSE}.}

\item{scale}{A boolean value specifying whether the two SRVFs should be
projected onto the Hilbert sphere before computing the distance. Defaults
to \code{FALSE}.}

\item{lambda}{A numeric value specifying the regularization parameter for the
optimal alignment. Defaults to \code{0}.}

\item{M}{An integer value specifying the number of points to use when
searching for the optimal rotation and alignment. Defaults to \code{200L}.}
}
\value{
A list with the following components:
\itemize{
\item \code{amplitude_distance}: A numeric value storing the amplitude distance
between the two SRVFs.
\item \code{phase_distance}: A numeric value storing the phase distance between the
two SRVFs.
\item \code{optimal_warping}: A function that takes a numeric vector \eqn{s} of values
in \eqn{[0, 1]} as input and returns the optimal warping function evaluated
at \eqn{s}.
\item \code{q1fun_modified}: A function that takes a numeric vector \eqn{s} of values
in \eqn{[0, 1]} as input and returns the first SRVF modified according to the
optimal alignment, rotation, and scaling.
\item \code{q2fun_modified}: A function that takes a numeric vector \eqn{s} of values
in \eqn{[0, 1]} as input and returns the second SRVF modified according to
the optimal alignment, rotation, and scaling.
}
}
\description{
Computes the distance between two shapes
}
\examples{
q1 <- curve2srvf(beta[, , 1, 1])
q2 <- curve2srvf(beta[, , 1, 2])
get_shape_distance(q1, q2)
}
