\name{horizFPCA}
\alias{horizFPCA}
\title{Horizontal Functional Principal Component Analysis}
\usage{
  horizFPCA(vec, mu, no, showplot = TRUE)
}
\arguments{
  \item{vec}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} of
  shooting vectors calculated using \code{\link{SqrtMean}}}

  \item{mu}{mean of warping functions}

  \item{no}{number of prinicpal components to extract}

  \item{showplot}{show plots of prinipal directions
  (default = T)}
}
\value{
  Returns a list containing \item{gam_pca}{warping
  functions principal directions} \item{psi_pca}{srvf
  principal directions} \item{latent}{latent values}
  \item{U}{eigenvectors}
}
\description{
  This function calculates vertical functional principal
  component analysis on aligned data
}
\examples{
data("simu_data")
gam = time_warping(f,time)$gam
out = SqrtMean(gam)
hfpca = horizFPCA(out$vec,out$mu,no = 3)
}
\references{
  Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron,
  J. S., May 2011. Registration of functional data using
  fisher-rao metric, arXiv:1103.3817v2 [math.ST].
}
\keyword{alignment}
\keyword{srvf}

