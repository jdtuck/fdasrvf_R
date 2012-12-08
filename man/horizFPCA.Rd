\name{horizFPCA}
\alias{horizFPCA}
\title{Horizontal Functional Principal Component Analysis}
\usage{
  horizFPCA(gam, no, showplot = TRUE)
}
\arguments{
  \item{gam}{matrix (\eqn{M} x \eqn{N}) of \eqn{M} of
  warping functions}

  \item{no}{number of prinicpal components to extract}

  \item{showplot}{show plots of prinipal directions
  (default = T)}
}
\value{
  Returns a list containing \item{gam_pca}{warping
  functions principal directions} \item{psi_pca}{srvf
  principal directions} \item{latent}{latent values}
  \item{U}{eigenvectors} \item{vec}{shooting vectors}
  \item{mu}{Karcher Mean}
}
\description{
  This function calculates vertical functional principal
  component analysis on aligned data
}
\examples{
data("simu_warp")
hfpca = horizFPCA(simu_warp$gam,no = 3)
}
\references{
  Tucker, J. D., Wu, W., Srivastava, A., Dec 2012.
  Generative Models for Function Data using Phase and
  Amplitude Separation, accepted to Computational
  Statistics and Data Analysis.
}
\keyword{alignment}
\keyword{srvf}

